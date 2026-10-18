#include "volume/lra_volume.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>
#include <lp_lib.h>

#include <random>

// GMP arbitrary-precision floats for the sampling / union bookkeeping. The
// default (--gmp) path stores points and evaluates membership tests in mpf at a
// arithmetic precision derived from lattice precision b; --fullgmp additionally runs
// the selected walk itself in this type (see gmp_walks.hpp).
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "logger.hpp"
#include "sampling/sphere.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/gmp_walks.hpp"

// cddlib (double precision) for polytope canonicalization, mirroring the
// Python prototype's use of pycddlib (src/polytope_operations.py::canonicalize).
#define set_int cdd_set_int
extern "C" {
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>
}
#undef set_int

namespace ttc
{
namespace
{

using TermIndexMap = std::unordered_map<cvc5::Term, std::size_t>;

// MpFloat (GMP float for the sampling / union bookkeeping) is defined in
// gmp_walks.hpp; its runtime precision is set globally via
// MpFloat::default_precision(digits) before sampling starts.

constexpr double kZeroTolerance = 1e-9;

// The intersection of one-sided coordinate bounds is a box. This deliberately
// rejects coupled facets, unbounded boxes, and degenerate boxes; it is an
// opt-in benchmark oracle, not a general polytope-volume algorithm.
double exactAxisAlignedBoxVolume(const Eigen::MatrixXd& A,
                                 const Eigen::VectorXd& b)
{
  const Eigen::Index dim = A.cols();
  std::vector<double> lower(dim, -std::numeric_limits<double>::infinity());
  std::vector<double> upper(dim, std::numeric_limits<double>::infinity());
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    Eigen::Index axis = -1;
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      if (std::abs(A(i, j)) <= 1e-12) continue;
      if (axis >= 0)
        throw std::runtime_error("--exact-box-volumes requires axis-aligned box facets");
      axis = j;
    }
    if (axis < 0)
    {
      if (b(i) < 0.0)
        throw std::runtime_error("--exact-box-volumes found an infeasible box");
      continue;
    }
    const double bound = b(i) / A(i, axis);
    if (!std::isfinite(bound))
      throw std::runtime_error("--exact-box-volumes found a nonfinite bound");
    if (A(i, axis) > 0)
      upper[axis] = std::min(upper[axis], bound);
    else
      lower[axis] = std::max(lower[axis], bound);
  }
  double volume = 1.0;
  for (Eigen::Index j = 0; j < dim; ++j)
  {
    if (!std::isfinite(lower[j]) || !std::isfinite(upper[j]) ||
        lower[j] >= upper[j])
      throw std::runtime_error("--exact-box-volumes requires a bounded, positive box");
    volume *= upper[j] - lower[j];
  }
  if (!std::isfinite(volume) || volume <= 0.0)
    throw std::runtime_error("--exact-box-volumes exceeds finite double volume range");
  return volume;
}

// Write the polytope {x : A x <= b} to `path` in the cdd / Avis-Fukuda
// H-representation format (readable by vinci, lrs, cddlib). Each cdd row is
// [b_i, -a_i0, -a_i1, ...], encoding b_i - a_i . x >= 0, i.e. a_i . x <= b_i.
void writePolytopeIne(const std::string& path, const Eigen::MatrixXd& A,
                      const Eigen::VectorXd& b)
{
  std::ofstream out(path);
  if (!out)
  {
    return;
  }
  const Eigen::Index rows = A.rows();
  const Eigen::Index cols = A.cols();
  out << "ttc_polytope\n";
  out << "H-representation\n";
  out << "begin\n";
  out << " " << rows << " " << (cols + 1) << " real\n";
  out << std::setprecision(17);
  for (Eigen::Index i = 0; i < rows; ++i)
  {
    out << " " << b(i);
    for (Eigen::Index j = 0; j < cols; ++j)
    {
      out << " " << -A(i, j);
    }
    out << "\n";
  }
  out << "end\n";
}

double parseRationalString(const std::string& value)
{
  if (value.empty())
  {
    return 0.0;
  }
  auto slash = value.find('/');
  if (slash == std::string::npos)
  {
    return std::stod(value);
  }
  double numerator = std::stod(value.substr(0, slash));
  double denominator = std::stod(value.substr(slash + 1));
  if (denominator == 0.0)
  {
    throw std::runtime_error("encountered zero denominator in rational value");
  }
  return numerator / denominator;
}

std::optional<double> getConstantValue(const cvc5::Term& term)
{
  if (term.isRealValue())
  {
    return parseRationalString(term.getRealValue());
  }
  if (term.isIntegerValue())
  {
    return parseRationalString(term.getIntegerValue());
  }
  switch (term.getKind())
  {
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return getConstantValue(term[0]);
    case cvc5::Kind::NEG:
      if (auto inner = getConstantValue(term[0]))
      {
        return -(*inner);
      }
      return std::nullopt;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        if (auto inner = getConstantValue(term[0]))
        {
          return -(*inner);
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto numerator = getConstantValue(term[0]);
        auto denominator = getConstantValue(term[1]);
        if (numerator && denominator && *denominator != 0.0)
        {
          return *numerator / *denominator;
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::MULT:
    {
      double product = 1.0;
      bool hasFactor = false;
      for (const auto& child : term)
      {
        auto factor = getConstantValue(child);
        if (!factor)
        {
          return std::nullopt;
        }
        product *= *factor;
        hasFactor = true;
      }
      if (hasFactor)
      {
        return product;
      }
      return std::nullopt;
    }
    default:
      break;
  }
  return std::nullopt;
}

bool accumulateLinear(const cvc5::Term& term,
                      double multiplier,
                      const TermIndexMap& index,
                      Eigen::VectorXd& coefficients,
                      double& constant)
{
  if (multiplier == 0.0)
  {
    return true;
  }

  if (term.isRealValue())
  {
    constant += multiplier * parseRationalString(term.getRealValue());
    return true;
  }
  if (term.isIntegerValue())
  {
    constant += multiplier * parseRationalString(term.getIntegerValue());
    return true;
  }

  if (term.getNumChildren() == 0)
  {
    if (term.getSort().isReal() || term.getSort().isInteger())
    {
      auto it = index.find(term);
      if (it == index.end())
      {
        Log(2) << "Unknown variable in linear term: " << term << std::endl;
        return false;
      }
      coefficients[it->second] += multiplier;
      return true;
    }
    return false;
  }

  switch (term.getKind())
  {
    case cvc5::Kind::ADD:
      for (const auto& child : term)
      {
        if (!accumulateLinear(child, multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
      }
      if (!accumulateLinear(term[0], multiplier, index, coefficients, constant))
      {
        return false;
      }
      for (std::size_t i = 1, n = term.getNumChildren(); i < n; ++i)
      {
        if (!accumulateLinear(term[i], -multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::NEG:
      return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return accumulateLinear(term[0], multiplier, index, coefficients, constant);
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto denom = getConstantValue(term[1]);
        if (!denom || *denom == 0.0)
        {
          return false;
        }
        return accumulateLinear(term[0], multiplier / *denom, index, coefficients,
                                constant);
      }
      return false;
    case cvc5::Kind::MULT:
    {
      double factor = 1.0;
      std::optional<cvc5::Term> expr;
      for (const auto& child : term)
      {
        auto constantValue = getConstantValue(child);
        if (constantValue)
        {
          factor *= *constantValue;
          continue;
        }
        if (expr.has_value())
        {
          return false;
        }
        expr = child;
      }
      if (!expr.has_value())
      {
        constant += multiplier * factor;
        return true;
      }
      return accumulateLinear(expr.value(), multiplier * factor, index,
                              coefficients, constant);
    }
    default:
      break;
  }
  return false;
}

struct Inequality
{
  Eigen::VectorXd coefficients;
  double bound = 0.0;  // Ax <= bound
};

struct PolytopeMatrices
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  bool empty = false;
};

bool appendLessEqualConstraint(const cvc5::Term& lhs,
                               const cvc5::Term& rhs,
                               const TermIndexMap& index,
                               std::size_t dimension,
                               std::vector<Inequality>& inequalities,
                               bool& empty)
{
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(static_cast<long>(dimension));
  double constant = 0.0;
  if (!accumulateLinear(lhs, 1.0, index, coeffs, constant) ||
      !accumulateLinear(rhs, -1.0, index, coeffs, constant))
  {
    return false;
  }
  double norm = coeffs.norm();
  if (norm <= kZeroTolerance)
  {
    if (constant <= kZeroTolerance)
    {
      return true;
    }
    empty = true;
    return true;
  }
  inequalities.push_back({std::move(coeffs), -constant});
  return true;
}

bool appendComparisonConstraints(const cvc5::Term& atom,
                                 bool positive,
                                 const TermIndexMap& index,
                                 std::size_t dimension,
                                 std::vector<Inequality>& inequalities,
                                 bool& empty)
{
  cvc5::Kind kind = atom.getKind();
  if (!positive)
  {
    switch (kind)
    {
      case cvc5::Kind::LEQ:
        kind = cvc5::Kind::GEQ;
        break;
      case cvc5::Kind::LT:
        kind = cvc5::Kind::GT;
        break;
      case cvc5::Kind::GEQ:
        kind = cvc5::Kind::LEQ;
        break;
      case cvc5::Kind::GT:
        kind = cvc5::Kind::LT;
        break;
      case cvc5::Kind::EQUAL:
        throw std::runtime_error("negated equality is not supported for LRA volume computation");
      default:
        break;
    }
  }

  switch (kind)
  {
    case cvc5::Kind::LEQ:
    case cvc5::Kind::LT:
      return appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::GEQ:
    case cvc5::Kind::GT:
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::EQUAL:
      if (!appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                     inequalities, empty))
      {
        return false;
      }
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    default:
      break;
  }
  return false;
}

bool collectConstraints(const cvc5::Term& term,
                        bool positive,
                        const TermIndexMap& index,
                        std::size_t dimension,
                        std::vector<Inequality>& inequalities,
                        bool& empty)
{
  if (empty)
  {
    return true;
  }

  if (term.getKind() == cvc5::Kind::NOT)
  {
    return collectConstraints(term[0], !positive, index, dimension, inequalities,
                              empty);
  }

  if (term.isBooleanValue())
  {
    bool value = term.getBooleanValue();
    bool satisfied = positive ? value : !value;
    if (!satisfied)
    {
      empty = true;
    }
    return true;
  }

  if (term.getKind() == cvc5::Kind::AND)
  {
    for (const auto& child : term)
    {
      if (!collectConstraints(child, positive, index, dimension, inequalities,
                              empty))
      {
        return false;
      }
    }
    return true;
  }

  return appendComparisonConstraints(term, positive, index, dimension,
                                     inequalities, empty);
}

PolytopeMatrices buildPolytopeMatrices(const Polytope& poly,
                                       const TermIndexMap& index,
                                       std::size_t dimension)
{
  PolytopeMatrices matrices;
  std::vector<Inequality> inequalities;
  inequalities.reserve(poly.termLiterals.size());
  bool empty = false;

  for (const auto& literal : poly.termLiterals)
  {
    if (!literal.value.has_value())
    {
      continue;
    }
    if (!collectConstraints(literal.term, literal.value.value(), index, dimension,
                            inequalities, empty))
    {
      throw std::runtime_error("unsupported literal encountered during LRA volume conversion");
    }
    if (empty)
    {
      break;
    }
  }

  if (empty)
  {
    matrices.empty = true;
    matrices.A = Eigen::MatrixXd::Zero(0, static_cast<long>(dimension));
    matrices.b = Eigen::VectorXd::Zero(0);
    return matrices;
  }

  if (inequalities.empty())
  {
    throw std::runtime_error("polytope has no linear constraints after abstraction");
  }

  matrices.A = Eigen::MatrixXd(static_cast<long>(inequalities.size()),
                               static_cast<long>(dimension));
  matrices.b = Eigen::VectorXd(static_cast<long>(inequalities.size()));

  for (std::size_t i = 0; i < inequalities.size(); ++i)
  {
    matrices.A.row(static_cast<long>(i)) = inequalities[i].coefficients.transpose();
    matrices.b[static_cast<long>(i)] = inequalities[i].bound;
  }

  return matrices;
}

// Canonicalizes the polytope {x : A x <= b} with cddlib, mirroring the Python
// prototype (src/polytope_operations.py::canonicalize -> cdd.matrix_canonicalize).
// Redundant inequalities are removed in place.  Returns false when the polytope
// is empty or has implicit equalities (i.e. it is lower-dimensional, so its
// full-dimensional volume is zero) -- the same cases where the Python code
// returns -1 / skips the polytope.  On any cddlib error the matrices are left
// unchanged and true is returned (fall back to the raw representation).
bool canonicalizePolytope(Eigen::MatrixXd& A, Eigen::VectorXd& b)
{
  const long m = A.rows();
  const long n = A.cols();
  if (m == 0)
  {
    return true;
  }

  dd_set_global_constants();

  // cddlib H-representation: row (g_0, g_1, ..., g_n) means
  // g_0 + g_1 x_1 + ... + g_n x_n >= 0.  For A x <= b that is b - A x >= 0,
  // so the row is (b_i, -A_i1, ..., -A_in).
  dd_MatrixPtr M = dd_CreateMatrix(m, n + 1);
  M->representation = dd_Inequality;
  for (long i = 0; i < m; ++i)
  {
    dd_set_d(M->matrix[i][0], b(i));
    for (long j = 0; j < n; ++j)
    {
      dd_set_d(M->matrix[i][j + 1], -A(i, j));
    }
  }

  dd_rowset implLin = nullptr;
  dd_rowset redset = nullptr;
  dd_rowindex newpos = nullptr;
  dd_ErrorType err = dd_NoError;
  dd_MatrixCanonicalize(&M, &implLin, &redset, &newpos, &err);

  bool keep = true;
  bool replace = false;
  Eigen::MatrixXd newA;
  Eigen::VectorXd newB;

  if (err != dd_NoError || M == nullptr)
  {
    // Leave the raw matrices untouched.
    keep = true;
  }
  else if (set_card(M->linset) > 0)
  {
    // Implicit equalities => measure-zero in the full dimension.
    keep = false;
  }
  else
  {
    long rows = M->rowsize;
    newA = Eigen::MatrixXd(rows, n);
    newB = Eigen::VectorXd(rows);
    for (long i = 0; i < rows; ++i)
    {
      newB(i) = dd_get_d(M->matrix[i][0]);
      for (long j = 0; j < n; ++j)
      {
        newA(i, j) = -dd_get_d(M->matrix[i][j + 1]);
      }
    }
    replace = true;
  }

  if (implLin != nullptr) set_free(implLin);
  if (redset != nullptr) set_free(redset);
  if (newpos != nullptr) free(newpos);
  if (M != nullptr) dd_FreeMatrix(M);
  dd_free_global_constants();

  if (replace)
  {
    A = std::move(newA);
    b = std::move(newB);
  }
  return keep;
}

// Stores the running set of sample points "X" used by the union-volume
// (Karp-Luby style) estimator.  Mirrors the behaviour of the reference Python
// implementation in src/cube_processor_nondis.py.  Templated on the coordinate
// type: `double` for the --nogmp path, `MpFloat` for the GMP paths (the
// membership test A x <= b is then evaluated in GMP, avoiding the double
// cancellation that can misclassify points near a facet).
template <typename Coord>
class SampleStore
{
 public:
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using Mat = Eigen::Matrix<Coord, Eigen::Dynamic, Eigen::Dynamic>;

  explicit SampleStore(double tolerance) : d_tolerance(tolerance) {}

  void addSample(const Vec& sample)
  {
    d_samples.push_back(sample);
  }

  // Removes every stored point that lies inside the polytope {x : A x <= b}.
  // Equivalent to the Python:  X = [s for s in X if not poly.is_in(s)].  The
  // half-space matrices are double (built once); they are cast to the
  // coordinate type so the dot products run at the stored precision.
  std::size_t removeInside(const Eigen::MatrixXd& Ad, const Eigen::VectorXd& bd)
  {
    if (Ad.rows() == 0)
    {
      return 0;
    }
    Mat A = Ad.template cast<Coord>();
    Vec b = bd.template cast<Coord>();
    Coord tol(d_tolerance);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec& sample) {
                               Vec diff = A * sample - b;
                               return (diff.array() <= tol).all();
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  // Keeps each stored point independently with probability 1/2.  Mirrors the
  // Python:  X = [s for s in X if random.random() >= 0.5].
  template <typename URNG>
  std::size_t downsampleHalf(URNG& rng)
  {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto oldSize = d_samples.size();
    auto it = std::remove_if(d_samples.begin(), d_samples.end(),
                             [&](const Vec&) {
                               return unit(rng) < 0.5;
                             });
    d_samples.erase(it, d_samples.end());
    return oldSize - d_samples.size();
  }

  std::size_t size() const
  {
    return d_samples.size();
  }

 private:
  double d_tolerance;
  std::vector<Vec> d_samples;
};

// GetPrecision in the revised manuscript: eta = epsilon/96. Work in logs
// so gamma/minRadius need not be representable as a double.
int precisionFromCubes(std::size_t dimension, std::size_t facetCount,
                       std::size_t numCubes, double minRadius, double epsilon)
{
  const double eta = epsilon / 96.0;
  const double logGamma = std::log(64.0) + std::log(double(numCubes))
      + std::log(double(dimension)) - std::log(eta)
      + 0.5 * std::log(std::log(4.0) + std::log(double(facetCount)) - std::log(eta));
  const double prec = std::max(0.0, std::ceil((logGamma - std::log(minRadius))
                                           / std::log(10.0)));
  if (!std::isfinite(prec) || prec > std::numeric_limits<int>::max() - 32)
  {
    throw std::runtime_error("lattice precision exceeds supported range");
  }
  return static_cast<int>(prec);
}

// Coordinate extrema bound P by a box, hence by a ball about center.
// Like the existing inner-ball LP, these bounds use practical double arithmetic.
double outerRadius(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& center)
{
  const int dim = static_cast<int>(A.cols());
  std::unique_ptr<lprec, decltype(&delete_lp)> lp(make_lp(0, dim), delete_lp);
  if (!lp) throw std::runtime_error("outer-radius LP allocation failed");
  set_verbose(lp.get(), CRITICAL);
  for (int j = 1; j <= dim; ++j) set_unbounded(lp.get(), j);
  std::vector<double> row(dim + 1, 0.0);
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    for (int j = 0; j < dim; ++j) row[j + 1] = A(i, j);
    if (!add_constraint(lp.get(), row.data(), LE, b(i)))
      throw std::runtime_error("outer-radius LP constraint failed");
  }
  double radius = 0.0;
  for (int j = 0; j < dim; ++j)
  {
    std::fill(row.begin(), row.end(), 0.0); row[j + 1] = 1.0;
    if (!set_obj_fn(lp.get(), row.data()))
      throw std::runtime_error("outer-radius LP objective failed");
    set_maxim(lp.get());
    if (solve(lp.get()) != OPTIMAL)
      throw std::runtime_error("polytope is unbounded or outer-radius LP failed");
    const double hi = get_objective(lp.get());
    set_minim(lp.get());
    if (solve(lp.get()) != OPTIMAL)
      throw std::runtime_error("polytope is unbounded or outer-radius LP failed");
    const double lo = get_objective(lp.get());
    radius = std::hypot(radius, std::max(std::abs(lo - center(j)),
                                        std::abs(hi - center(j))));
  }
  if (!std::isfinite(radius) || radius <= 0.0)
    throw std::runtime_error("invalid polytope outer radius");
  return radius;
}

// GenerateSamples: inflate facets, independently restart the selected walk,
// stochastically round to 10^-b Z^n, and reject against the original body.
// Walk length and Volesti convergence remain practical, uncertified oracles.
template <typename Coord, typename GenFn>
std::vector<Eigen::Matrix<Coord, Eigen::Dynamic, 1>> sampleLattice(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
    int precision, double epsilon, double delta, unsigned seed,
    GenFn&& generate, double& samplingTime)
{
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using std::floor;
  using std::pow;
  const Coord scale = pow(Coord(10), precision);
  const double spacing = std::pow(10.0, -precision);
  if (spacing == 0.0)
    throw std::runtime_error("lattice spacing is below double geometry resolution");
  const double kappa = std::sqrt(2.0 * std::log(16.0 / epsilon))
                       + std::sqrt(2.0 * std::log(double(A.rows())));
  Eigen::VectorXd inflated = b;
  for (Eigen::Index i = 0; i < A.rows(); ++i)
    inflated(i) += spacing * kappa * A.row(i).norm();
  const auto Am = A.template cast<Coord>().eval();
  const Vec bm = b.template cast<Coord>();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::vector<Vec> out;
  out.reserve(static_cast<std::size_t>(n));
  // Under the paper's oracle assumptions, rejection probability <= epsilon/3.
  // Budget a failed sampling batch within delta; report failure, never fewer points.
  const int maxBatches = std::max(1, static_cast<int>(std::ceil(
      (std::log(delta) - std::log(double(std::max(1LL, n)))) / std::log(epsilon / 3.0))));
  const double begun = Log.elapsed();
  for (int batch = 0; out.size() < static_cast<std::size_t>(n); ++batch)
  {
    if (batch >= maxBatches)
      throw std::runtime_error("lattice rejection sampler exceeded its retry budget");
    double ignoredTime = 0.0;
    auto points = generate(A, inflated, n - static_cast<long long>(out.size()), ignoredTime);
    for (const auto& point : points)
    {
      Vec rounded(point.size());
      for (Eigen::Index j = 0; j < point.size(); ++j)
      {
        const Coord raw = Coord(point(j)) * scale;
        Coord base = floor(raw);
        const Coord fraction = raw - base;
        if (unit(rng) < static_cast<double>(fraction)) base += 1;
        rounded(j) = base / scale;
      }
      if (((Am * rounded - bm).array() <= Coord(0)).all())
        out.push_back(std::move(rounded));
    }
  }
  samplingTime = Log.elapsed() - begun;
  return out;
}

// A polytope kept after phase 1, with its (double) half-space matrices and the
// volesti volume estimate.
struct PolytopeVolume
{
  Eigen::MatrixXd A;
  Eigen::VectorXd b;
  double volume = 0.0;
};

// Phase 2: streaming union-volume (Karp-Luby style) estimation, templated on the
// coordinate type used to store the surviving sample set `X`.  `p` is the
// running sampling probability; the union volume is |X| / (p*C), C=T*10^(b*n).
// `generate(A, b, n, samplingTime)` draws `n` fresh points inside {x : A x <= b}
// (the only step that differs between the double and GMP walk backends).
template <typename Coord, typename GenFn>
void runUnionAlgorithm(const std::vector<PolytopeVolume>& kept, double thresh,
                       unsigned seed, double log2Scale, GenFn&& generate,
                       VolumeComputationResult& result,
                       const std::function<void(const VolumeComputationRow&)>& onRow)
{
  std::mt19937 auxRng(seed);  // drives Poisson draws and down-sampling
  SampleStore<Coord> store(0.0);

  // Track log2(p*C), avoiding overflow of C and underflow of p.
  double log2Density = log2Scale;
  std::size_t rowIndex = 0;
  std::size_t totalGenerated = 0;
  std::size_t totalDeleted = 0;
  double totalSamplingTime = 0.0;

  for (const auto& entry : kept)
  {
    const double volume = entry.volume;

    // Drop stored points that lie inside the current polytope.
    std::size_t deleted = store.removeInside(entry.A, entry.b);

    // Shrink the sampling probability until p*C*volume < thresh.
    while (log2Density + std::log2(volume) >= std::log2(thresh))
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
    }

    // Draw the number of fresh samples; keep |X| + N within the threshold.
    std::poisson_distribution<long long> poisson(std::exp2(log2Density + std::log2(volume)));
    long long n = poisson(auxRng);
    while (static_cast<double>(n) + static_cast<double>(store.size()) > thresh)
    {
      store.downsampleHalf(auxRng);
      log2Density -= 1.0;
      n = std::binomial_distribution<long long>(n, 0.5)(auxRng);
    }

    std::size_t sampleCount = 0;
    double samplingTime = 0.0;
    if (n > 0)
    {
      auto samples = generate(entry.A, entry.b, n, samplingTime);
      for (auto& sample : samples)
      {
        store.addSample(sample);
      }
      sampleCount = samples.size();
    }

    totalGenerated += sampleCount;
    totalDeleted += deleted;
    totalSamplingTime += samplingTime;
    VolumeComputationRow row;
    row.index = ++rowIndex;
    row.volume = volume;
    row.samplesDeleted = deleted;
    row.totalSamples = store.size();
    row.samplesGenerated = sampleCount;
    result.rows.push_back(row);
    if (onRow)
    {
      onRow(row);
    }
  }

  result.finalSampleCount = store.size();
  result.volumeEstimate =
      store.size() ? std::exp2(std::log2(double(store.size())) - log2Density) : 0.0;
  result.totalSamplesGenerated = totalGenerated;
  result.totalSamplesDeleted = totalDeleted;
  result.samplingTime = totalSamplingTime;
}

std::vector<cvc5::Term> deduplicateRealVariables(
    const std::vector<cvc5::Term>& variables)
{
  std::vector<cvc5::Term> unique;
  unique.reserve(variables.size());
  std::unordered_set<cvc5::Term> seen;
  for (const auto& term : variables)
  {
    if (!term.getSort().isReal())
    {
      continue;
    }
    if (seen.insert(term).second)
    {
      unique.push_back(term);
    }
  }
  return unique;
}

}  // namespace

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const VolumeOptions& options,
    const std::function<void(const VolumeComputationRow&)>& onRow)
{
  VolumeComputationResult result;
  if (!std::isfinite(options.epsilon) || options.epsilon <= 0.0 || options.epsilon > 1.0
      || !std::isfinite(options.delta) || options.delta <= 0.0 || options.delta >= 1.0)
    throw std::invalid_argument("LRA volume requires 0 < epsilon <= 1 and 0 < delta < 1");
  if (options.sampleWalkLength == 0 || options.volumeWalkLength == 0 || options.precision < 0)
    throw std::invalid_argument("volume walk lengths and arithmetic precision must be positive or automatic");
  auto uniqueVars = deduplicateRealVariables(realVariables);
  if (uniqueVars.empty())
  {
    throw std::runtime_error("LRA volume computation requires at least one real variable");
  }

  std::size_t dimension = uniqueVars.size();
  TermIndexMap index;
  index.reserve(dimension);
  for (std::size_t i = 0; i < uniqueVars.size(); ++i)
  {
    index.emplace(uniqueVars[i], i);
  }

  using Kernel = Cartesian<double>;
  using Point = Kernel::Point;
  using HPolytope = HPolytope<Point>;
  // Volume RNG matches the volesti `volume` tool: a fixed compile-time seed (3),
  // so the cooling-gaussians estimate is deterministic.
  using VolumeRNG = BoostRandomNumberGenerator<boost::mt19937, double, 3>;
  // Sampling RNG matches the volesti `sample` tool: explicitly seeded.
  using SampleRNG = BoostRandomNumberGenerator<boost::mt19937, double>;
  // Candidate sampling walks, selected at run time by options.samplerWalk.
  using BilliardW = BilliardWalk::template Walk<HPolytope, SampleRNG>;
  using AccelBilliardW =
      AcceleratedBilliardWalk::template Walk<HPolytope, SampleRNG>;
  using BallW = BallWalk::template Walk<HPolytope, SampleRNG>;
  using RdhrW = RDHRWalk::template Walk<HPolytope, SampleRNG>;
  using CdhrW = CDHRWalk::template Walk<HPolytope, SampleRNG>;

  // Algorithm 1 in the revised manuscript.
  const double kEpsilon = options.epsilon;
  const double kDelta = options.delta;
  const unsigned kSeed = options.seed;
  // Per-polytope sampling streams.  `generate` is invoked once per polytope by
  // runUnionAlgorithm, and each invocation seeds its walk from scratch, so a
  // single base seed would replay the identical direction sequence for every
  // polytope and correlate the sample sets the union estimator combines.
  // Derive a distinct stream per call instead, mixing with splitmix64 so
  // consecutive polytopes do not get near-adjacent mt19937 states.
  unsigned streamIndex = 0;
  auto streamSeed = [&kSeed, &streamIndex]() {
    std::uint64_t z = kSeed + 0x9E3779B97F4A7C15ULL * (++streamIndex);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return static_cast<unsigned>(z ^ (z >> 31));
  };
  const double mvcEps = kEpsilon / 12.0;
  const double numCubes = static_cast<double>(polytopes.size());
  const double oracleDelta = kDelta / (4.0 * std::max(1.0, numCubes));
  const double thresh = std::max(
      24.0 * (std::log(24.0) - std::log(kDelta)) / ((1.0 - mvcEps) * mvcEps * mvcEps),
      15.0 * (std::log(6.0) - std::log(kDelta) + std::log(std::max(1.0, numCubes))));
  if (!std::isfinite(thresh) || thresh >= double(std::numeric_limits<long long>::max() / 2)
      || oracleDelta <= 0.0)
    throw std::invalid_argument("epsilon/delta require an unsupported sample threshold");

  // The volesti `volume` tool uses walk length 10 + dim/10; the `sample` tool
  // uses a fixed walk length of 10. Both are overridable via --walklen-vol /
  // --walklen-samp (a negative option value keeps the volesti default).
  const unsigned volumeWalkLength =
      options.volumeWalkLength >= 0
          ? static_cast<unsigned>(options.volumeWalkLength)
          : static_cast<unsigned>(10 + dimension / 10);
  const unsigned sampleWalkLength =
      options.sampleWalkLength >= 0
          ? static_cast<unsigned>(options.sampleWalkLength)
          : 10;

  Log(2) << "Starting union algorithm, threshold: " << thresh << std::endl;
  // Volesti's practical cooling API accepts epsilon but no confidence delta.
  Log(2) << "Per-polytope oracle targets: epsilon=" << mvcEps
         << ", delta=" << oracleDelta << " (practical Volesti oracle)" << std::endl;

  // Phase 1: build the half-space matrices and compute every volume up front,
  // dropping polytopes whose (full-dimensional) volume is zero.
  std::vector<PolytopeVolume> kept;
  kept.reserve(polytopes.size());

  std::size_t numZeroVolume = 0;
  double totalVolumeTime = 0.0;
  std::size_t facetTotal = 0;
  double minRadius = std::numeric_limits<double>::infinity();
  std::size_t dumpIndex = 0;  // running counter for --dump-ine file names
  for (const auto& poly : polytopes)
  {
    PolytopeMatrices matrices = buildPolytopeMatrices(poly, index, dimension);
    if (matrices.empty)
    {
      ++numZeroVolume;
      continue;
    }
    if (!matrices.A.allFinite() || !matrices.b.allFinite())
      throw std::runtime_error("polytope coefficients exceed double geometry range");
    // Canonicalize (remove redundant constraints, detect equalities) with
    // cddlib, exactly as the Python prototype does before volume/sampling.
    // --no-cdd-simp skips this and feeds volesti the raw half-spaces.
    long rawFacets = matrices.A.rows();
    if (options.cddSimplify && !canonicalizePolytope(matrices.A, matrices.b))
    {
      ++numZeroVolume;
      continue;
    }
    if (!options.dumpInePrefix.empty())
    {
      std::string path = options.dumpInePrefix + "_cube" +
                         std::to_string(++dumpIndex) + ".ine";
      writePolytopeIne(path, matrices.A, matrices.b);
      Log(1) << "wrote polytope to " << path << std::endl;
    }
    HPolytope hpoly(static_cast<unsigned>(dimension), matrices.A, matrices.b);
    const auto innerBall = hpoly.ComputeInnerBall();
    if (innerBall.second <= 0.0)
    {
      ++numZeroVolume;
      continue;
    }
    if (!std::isfinite(innerBall.second))
      throw std::runtime_error("polytope has no finite inner radius");
    const double radius = outerRadius(matrices.A, matrices.b, innerBall.first.getCoefficients());
    double volumeComputeStart = Log.elapsed();
    double volume = options.exactBoxVolumes
        ? exactAxisAlignedBoxVolume(matrices.A, matrices.b)
        : volume_cooling_gaussians<GaussianBallWalk, VolumeRNG>(
              hpoly, mvcEps, volumeWalkLength);
    double volumeComputeElapsed = Log.elapsed() - volumeComputeStart;
    totalVolumeTime += volumeComputeElapsed;
    Log(2) << "polytope volume " << volume << " (" << matrices.A.rows()
           << " facets, was " << rawFacets << ", dim " << dimension << ") in "
           << volumeComputeElapsed << "s" << std::endl;
    if (!std::isfinite(volume) || volume <= 0.0)
      throw std::runtime_error("volume oracle failed on a positive-dimensional polytope");
    const double lower = (1.0 - mvcEps) * std::pow(2.0 * innerBall.second / std::sqrt(double(dimension)), double(dimension));
    const double upper = (1.0 + mvcEps) * std::pow(2.0 * radius, double(dimension));
    if (std::isnan(lower) || std::isnan(upper) || lower > upper)
      throw std::runtime_error("invalid polytope volume clipping bounds");
    if (!options.exactBoxVolumes)
      volume = std::clamp(volume, lower, upper);
    if (!std::isfinite(volume) || volume <= 0.0)
      throw std::runtime_error("clipped polytope volume exceeds double range");
    minRadius = std::min(minRadius, innerBall.second);
    facetTotal += static_cast<std::size_t>(matrices.A.rows());
    kept.push_back({std::move(matrices.A), std::move(matrices.b), volume});
  }

  if (numZeroVolume > 0)
  {
    Log(2) << "Skipping " << numZeroVolume << " polytopes out of "
           << polytopes.size() << " where volume is zero" << std::endl;
  }

  result.volumeComputationTime = totalVolumeTime;
  if (kept.empty()) return result;
  // Lattice spacing and GMP arithmetic precision are distinct quantities.
  const int precision = precisionFromCubes(dimension, facetTotal, polytopes.size(),
                                          minRadius, kEpsilon);
  result.latticePrecision = precision;
  const double log2Scale = std::log2(thresh) + double(precision) * double(dimension) * std::log2(10.0);
  if (options.gmpMode != GmpMode::None)
  {
    result.samplingPrecision = std::max(options.precision, precision + 32);
    MpFloat::default_precision(static_cast<unsigned>(result.samplingPrecision));
  }

  SampleRNG sampleRng(static_cast<int>(dimension));

  // Restart the selected walk at one fixed center before every proposal.
  auto generateDouble =
      [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
          double& samplingTime) -> std::vector<Eigen::VectorXd> {
    HPolytope hpoly(static_cast<unsigned>(dimension), A, b);
    const Point center = hpoly.ComputeInnerBall().first;
    PushBackWalkPolicy policy;
    std::list<Point> samples;
    sampleRng.set_seed(streamSeed());
    double samplingStart = Log.elapsed();
    for (long long i = 0; i < n; ++i)
    {
      Point start = center;
      switch (options.samplerWalk)
      {
        case SamplerWalk::AcceleratedBilliard:
          RandomPointGenerator<AccelBilliardW>::apply(
              hpoly, start, 1, sampleWalkLength, samples, policy, sampleRng);
          break;
        case SamplerWalk::Ball:
          RandomPointGenerator<BallW>::apply(
              hpoly, start, 1, sampleWalkLength, samples, policy, sampleRng);
          break;
        case SamplerWalk::RDHR:
          RandomPointGenerator<RdhrW>::apply(
              hpoly, start, 1, sampleWalkLength, samples, policy, sampleRng);
          break;
        case SamplerWalk::CDHR:
          RandomPointGenerator<CdhrW>::apply(
              hpoly, start, 1, sampleWalkLength, samples, policy, sampleRng);
          break;
        case SamplerWalk::Billiard:
        default:
          RandomPointGenerator<BilliardW>::apply(
              hpoly, start, 1, sampleWalkLength, samples, policy, sampleRng);
          break;
      }
    }
    samplingTime = Log.elapsed() - samplingStart;
    std::vector<Eigen::VectorXd> out;
    out.reserve(samples.size());
    for (const Point& sample : samples)
    {
      out.push_back(sample.getCoefficients());
    }
    return out;
  };

  if (options.gmpMode == GmpMode::Full &&
      options.samplerWalk == SamplerWalk::AcceleratedBilliard)
  {
    Log(0) << "WARNING: --fullgmp has no accelerated-billiard implementation; "
              "using the GMP billiard walk instead"
           << std::endl;
  }

  switch (options.gmpMode)
  {
    case GmpMode::None:
    {
      auto generate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                          long long n, double& time) {
        return sampleLattice<double>(A, b, n, precision, mvcEps, oracleDelta,
                                     streamSeed(), generateDouble, time);
      };
      runUnionAlgorithm<double>(kept, thresh, kSeed, log2Scale, generate, result,
                                onRow);
      break;
    }
    case GmpMode::PointRepr:
    {
      // Double proposals, with lattice rounding and membership in GMP.
      auto generate =
          [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
              double& samplingTime) {
            return sampleLattice<MpFloat>(A, b, n, precision, mvcEps, oracleDelta,
                                          streamSeed(), generateDouble, samplingTime);
          };
      runUnionAlgorithm<MpFloat>(kept, thresh, kSeed, log2Scale, generate, result, onRow);
      break;
    }
    case GmpMode::Full:
    {
      // Run the walk itself in GMP.  The Chebyshev centre and inner radius
      // (used to size the walk step) come from the double inner ball.  Billiard,
      // ball, RDHR and CDHR have GMP implementations; accelerated-billiard was
      // warned about above and falls back to the GMP billiard walk.
      const SamplerWalk walk = options.samplerWalk;
      auto generate =
          [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
              double& samplingTime) {
            HPolytope hpoly(static_cast<unsigned>(dimension), A, b);
            auto ball = hpoly.ComputeInnerBall();
            Eigen::VectorXd center = ball.first.getCoefficients();
            double radius = ball.second;
            double samplingStart = Log.elapsed();
            // One stream advance per polytope, drawn before the switch so the
            // sequence of seeds does not depend on which walk is selected.
            const unsigned walkSeed = streamSeed();
            std::vector<MpVector> pts;
            switch (walk)
            {
              case SamplerWalk::Ball:
                pts = sampleGmpBallWalk(A, b, center, radius, n,
                                        sampleWalkLength, walkSeed, true);
                break;
              case SamplerWalk::RDHR:
                pts = sampleGmpRDHR(A, b, center, radius, n, sampleWalkLength,
                                    walkSeed, true);
                break;
              case SamplerWalk::CDHR:
                pts = sampleGmpCDHR(A, b, center, radius, n, sampleWalkLength,
                                    walkSeed, true);
                break;
              case SamplerWalk::AcceleratedBilliard:
              case SamplerWalk::Billiard:
              default:
                pts = sampleGmpBilliard(A, b, center, radius, n,
                                        sampleWalkLength, walkSeed, true);
                break;
            }
            samplingTime = Log.elapsed() - samplingStart;
            return pts;
          };
      auto latticeGenerate = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                 long long n, double& time) {
        return sampleLattice<MpFloat>(A, b, n, precision, mvcEps, oracleDelta,
                                      streamSeed(), generate, time);
      };
      runUnionAlgorithm<MpFloat>(kept, thresh, kSeed, log2Scale, latticeGenerate, result, onRow);
      break;
    }
  }

  result.volumeComputationTime = totalVolumeTime;
  return result;
}

}  // namespace ttc

