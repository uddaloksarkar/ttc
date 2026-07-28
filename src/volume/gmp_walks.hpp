#pragma once

// Self-contained random walks (billiard, ball, and hit-and-run) in GMP (boost
// mpf_float) arithmetic, used by the LRA volume engine's --fullgmp mode.
// volesti's templated walks cannot be instantiated on mpf_float without an
// extensive port (boost.random does not support variable-precision floats, and
// several volesti routines narrow NT to int/double), so we implement the walks
// directly here.  The geometry is the same one volesti uses (billiard
// reflection / ball-walk accept-reject / hit-and-run chord sampling); only the
// arithmetic differs -- every dot product, reflection, chord intersection and
// membership test runs at the configured mpf precision, so points near a facet
// are placed without double cancellation.

#include <cstddef>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace ttc
{

// GMP float for the sampling / union bookkeeping.  Runtime precision is set
// globally via MpFloat::default_precision(digits) before sampling starts.
using MpFloat = boost::multiprecision::mpf_float;
using MpVector = Eigen::Matrix<MpFloat, Eigen::Dynamic, 1>;

// Draws `n` points from the billiard walk inside {x : A x <= b}, each `n`
// separated by `walkLength` reflection steps, starting from `startCenter` (the
// double Chebyshev centre, converted to mpf).  `innerRadius` sizes the walk
// step.  Direction draws use a double Mersenne twister seeded with `seed`
// (randomness need not exceed double precision); all geometry is in mpf.
inline std::vector<MpVector> sampleGmpBilliard(const Eigen::MatrixXd& Ad,
                                               const Eigen::VectorXd& bd,
                                               const Eigen::VectorXd& startCenter,
                                               double innerRadius, long long n,
                                               unsigned walkLength, unsigned seed)
{
  const Eigen::Index m = Ad.rows();
  const Eigen::Index dim = Ad.cols();
  std::vector<MpVector> out;
  if (n <= 0 || m == 0 || dim == 0)
  {
    return out;
  }

  Eigen::Matrix<MpFloat, Eigen::Dynamic, Eigen::Dynamic> Am = Ad.cast<MpFloat>();
  MpVector bm = bd.cast<MpFloat>();
  // Per-row squared norm <a_i, a_i> for the reflection, precomputed once.
  MpVector rowSqNorm(m);
  for (Eigen::Index i = 0; i < m; ++i)
  {
    rowSqNorm(i) = Am.row(i).dot(Am.row(i));
  }

  // Walk step length, mirroring volesti's diameter estimate 2*sqrt(dim)*r.
  MpFloat len = MpFloat(2) * sqrt(MpFloat(static_cast<long>(dim)))
                * MpFloat(innerRadius);
  if (len <= 0)
  {
    len = MpFloat(1);
  }
  const MpFloat dl(0.995);  // back off slightly from each facet (volesti)
  const MpFloat eps("1e-12");

  std::mt19937 rng(seed);
  std::normal_distribution<double> ndist(0.0, 1.0);
  std::uniform_real_distribution<double> udist(0.0, 1.0);

  MpVector p = startCenter.cast<MpFloat>();
  MpVector v(dim);
  out.reserve(static_cast<std::size_t>(n));

  auto sampleDirection = [&]() {
    MpFloat sq(0);
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) = MpFloat(ndist(rng));
      sq += v(j) * v(j);
    }
    MpFloat norm = sqrt(sq);
    if (norm <= 0)
    {
      norm = MpFloat(1);
    }
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) /= norm;
    }
  };

  const long maxReflections = 50 * static_cast<long>(dim);
  for (long long pt = 0; pt < n; ++pt)
  {
    for (unsigned step = 0; step < walkLength; ++step)
    {
      MpFloat budget = MpFloat(udist(rng)) * len;
      sampleDirection();
      MpVector p0 = p;
      long it = 0;
      for (; it < maxReflections; ++it)
      {
        // First positive facet along v: t_i = (b_i - a_i.p) / (a_i.v).
        MpVector av = Am * v;
        MpVector ap = Am * p;
        MpFloat minPlus(0);
        Eigen::Index facet = -1;
        for (Eigen::Index i = 0; i < m; ++i)
        {
          if (av(i) > eps)
          {
            MpFloat t = (bm(i) - ap(i)) / av(i);
            if (facet < 0 || t < minPlus)
            {
              minPlus = t;
              facet = i;
            }
          }
        }
        if (facet < 0)
        {
          break;  // unbounded along v (should not happen for a bounded body)
        }
        if (budget <= minPlus)
        {
          p += budget * v;
          break;
        }
        MpFloat advance = dl * minPlus;
        p += advance * v;
        budget -= advance;
        // Reflect v about the facet normal a: v -= 2 <v,a>/<a,a> a.
        MpFloat coeff = MpFloat(2) * Am.row(facet).dot(v) / rowSqNorm(facet);
        v -= coeff * Am.row(facet).transpose();
      }
      if (it == maxReflections)
      {
        p = p0;  // numerical trouble: discard this step
      }
    }
    out.push_back(p);
  }
  return out;
}

// Draws `n` points from the ball walk inside {x : A x <= b}, mirroring volesti's
// uniform BallWalk (uniform_ball_walk.hpp) but in mpf arithmetic.  Each step
// proposes y = p + (uniform point in the d-ball of radius delta) and accepts it
// iff A y <= b (metropolis for the uniform target).  The step radius matches
// volesti's compute_delta: delta = 4 * innerRadius / sqrt(dim).  Directions and
// the radial fraction are drawn in double (randomness need not exceed double
// precision); the proposal and the membership test A y <= b run in mpf, so a
// point is accepted/rejected without double cancellation near a facet.
inline std::vector<MpVector> sampleGmpBallWalk(const Eigen::MatrixXd& Ad,
                                               const Eigen::VectorXd& bd,
                                               const Eigen::VectorXd& startCenter,
                                               double innerRadius, long long n,
                                               unsigned walkLength, unsigned seed)
{
  const Eigen::Index m = Ad.rows();
  const Eigen::Index dim = Ad.cols();
  std::vector<MpVector> out;
  if (n <= 0 || m == 0 || dim == 0)
  {
    return out;
  }

  Eigen::Matrix<MpFloat, Eigen::Dynamic, Eigen::Dynamic> Am = Ad.cast<MpFloat>();
  MpVector bm = bd.cast<MpFloat>();

  // Ball-walk step radius, mirroring volesti's compute_delta = 4 r / sqrt(dim).
  MpFloat delta = MpFloat(4) * MpFloat(innerRadius)
                  / sqrt(MpFloat(static_cast<long>(dim)));
  if (delta <= 0)
  {
    delta = MpFloat(1);
  }
  const double invDim = 1.0 / static_cast<double>(dim);

  std::mt19937 rng(seed);
  std::normal_distribution<double> ndist(0.0, 1.0);
  std::uniform_real_distribution<double> udist(0.0, 1.0);

  MpVector p = startCenter.cast<MpFloat>();
  MpVector v(dim);
  out.reserve(static_cast<std::size_t>(n));

  // Uniform unit direction (normalized Gaussian), same convention as the
  // billiard walk above.
  auto sampleDirection = [&]() {
    MpFloat sq(0);
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) = MpFloat(ndist(rng));
      sq += v(j) * v(j);
    }
    MpFloat norm = sqrt(sq);
    if (norm <= 0)
    {
      norm = MpFloat(1);
    }
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) /= norm;
    }
  };

  for (long long pt = 0; pt < n; ++pt)
  {
    for (unsigned step = 0; step < walkLength; ++step)
    {
      sampleDirection();
      // Uniform radius in the d-ball: r = delta * U^(1/dim), U ~ U(0,1)
      // (GetPointInDsphere in volesti/sampling/sphere.hpp).
      MpFloat radius = delta * MpFloat(std::pow(udist(rng), invDim));
      MpVector y = p + radius * v;
      // Accept iff y stays inside every half-space (A y <= b).
      MpVector ay = Am * y;
      bool inside = true;
      for (Eigen::Index i = 0; i < m; ++i)
      {
        if (ay(i) > bm(i))
        {
          inside = false;
          break;
        }
      }
      if (inside)
      {
        p = y;
      }
    }
    out.push_back(p);
  }
  return out;
}

// Chord of the line {p + t v : t in R} inside {x : A x <= b}, as [tMin, tMax].
// `av` = A v and `ap` = A p are precomputed.  Returns false if the chord is
// empty/unbounded (degenerate direction) so the caller can skip the step.
inline bool gmpLineChord(const MpVector& av, const MpVector& ap,
                         const MpVector& bm, const MpFloat& eps, MpFloat& tMin,
                         MpFloat& tMax)
{
  const Eigen::Index m = av.size();
  bool haveMin = false, haveMax = false;
  for (Eigen::Index i = 0; i < m; ++i)
  {
    // Row i requires t * av(i) <= bm(i) - ap(i).
    if (av(i) > eps)  // upper bound on t
    {
      MpFloat t = (bm(i) - ap(i)) / av(i);
      if (!haveMax || t < tMax) { tMax = t; haveMax = true; }
    }
    else if (av(i) < -eps)  // lower bound on t
    {
      MpFloat t = (bm(i) - ap(i)) / av(i);
      if (!haveMin || t > tMin) { tMin = t; haveMin = true; }
    }
  }
  return haveMin && haveMax && tMax > tMin;
}

// Draws `n` points from the random-directions hit-and-run walk (RDHRWalk in
// volesti) inside {x : A x <= b}, in mpf arithmetic.  Each step picks a uniform
// direction v, computes the chord [tMin, tMax] of {p + t v} within the polytope,
// and jumps to a uniform point p + t v, t ~ U(tMin, tMax).  `innerRadius` is
// unused (hit-and-run needs no step size); kept for a uniform call signature.
inline std::vector<MpVector> sampleGmpRDHR(const Eigen::MatrixXd& Ad,
                                           const Eigen::VectorXd& bd,
                                           const Eigen::VectorXd& startCenter,
                                           double /*innerRadius*/, long long n,
                                           unsigned walkLength, unsigned seed)
{
  const Eigen::Index m = Ad.rows();
  const Eigen::Index dim = Ad.cols();
  std::vector<MpVector> out;
  if (n <= 0 || m == 0 || dim == 0)
  {
    return out;
  }

  Eigen::Matrix<MpFloat, Eigen::Dynamic, Eigen::Dynamic> Am = Ad.cast<MpFloat>();
  MpVector bm = bd.cast<MpFloat>();
  const MpFloat eps("1e-12");

  std::mt19937 rng(seed);
  std::normal_distribution<double> ndist(0.0, 1.0);
  std::uniform_real_distribution<double> udist(0.0, 1.0);

  MpVector p = startCenter.cast<MpFloat>();
  MpVector v(dim);
  out.reserve(static_cast<std::size_t>(n));

  auto sampleDirection = [&]() {
    MpFloat sq(0);
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) = MpFloat(ndist(rng));
      sq += v(j) * v(j);
    }
    MpFloat norm = sqrt(sq);
    if (norm <= 0)
    {
      norm = MpFloat(1);
    }
    for (Eigen::Index j = 0; j < dim; ++j)
    {
      v(j) /= norm;
    }
  };

  for (long long pt = 0; pt < n; ++pt)
  {
    for (unsigned step = 0; step < walkLength; ++step)
    {
      sampleDirection();
      MpVector av = Am * v;
      MpVector ap = Am * p;
      MpFloat tMin(0), tMax(0);
      if (!gmpLineChord(av, ap, bm, eps, tMin, tMax))
      {
        continue;  // degenerate direction; keep p and try again
      }
      MpFloat lambda = tMin + MpFloat(udist(rng)) * (tMax - tMin);
      p += lambda * v;
    }
    out.push_back(p);
  }
  return out;
}

// Draws `n` points from the coordinate-directions hit-and-run walk (CDHRWalk in
// volesti) inside {x : A x <= b}, in mpf arithmetic.  Identical to RDHR except
// the direction is a uniformly chosen coordinate axis e_c, so the chord is
// computed from column c of A.  `innerRadius` is unused.
inline std::vector<MpVector> sampleGmpCDHR(const Eigen::MatrixXd& Ad,
                                           const Eigen::VectorXd& bd,
                                           const Eigen::VectorXd& startCenter,
                                           double /*innerRadius*/, long long n,
                                           unsigned walkLength, unsigned seed)
{
  const Eigen::Index m = Ad.rows();
  const Eigen::Index dim = Ad.cols();
  std::vector<MpVector> out;
  if (n <= 0 || m == 0 || dim == 0)
  {
    return out;
  }

  Eigen::Matrix<MpFloat, Eigen::Dynamic, Eigen::Dynamic> Am = Ad.cast<MpFloat>();
  MpVector bm = bd.cast<MpFloat>();
  const MpFloat eps("1e-12");

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> udist(0.0, 1.0);
  std::uniform_int_distribution<int> cdist(0, static_cast<int>(dim) - 1);

  MpVector p = startCenter.cast<MpFloat>();
  // ap = A p, maintained incrementally: only coordinate c changes per step.
  MpVector ap = Am * p;
  out.reserve(static_cast<std::size_t>(n));

  for (long long pt = 0; pt < n; ++pt)
  {
    for (unsigned step = 0; step < walkLength; ++step)
    {
      const Eigen::Index c = static_cast<Eigen::Index>(cdist(rng));
      // Chord along axis e_c: av is column c of A.
      const MpVector av = Am.col(c);
      MpFloat tMin(0), tMax(0);
      if (!gmpLineChord(av, ap, bm, eps, tMin, tMax))
      {
        continue;
      }
      MpFloat t = tMin + MpFloat(udist(rng)) * (tMax - tMin);
      p(c) += t;
      ap += t * av;  // keep A p in sync with the single-coordinate move
    }
    out.push_back(p);
  }
  return out;
}

}  // namespace ttc
