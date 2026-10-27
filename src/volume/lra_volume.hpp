#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include <cvc5/cvc5.h>

#include "volume/allsat_volume.hpp"

namespace ttc
{

struct VolumeComputationRow
{
  std::size_t index = 0;
  double volume = 0.0;
  std::size_t samplesDeleted = 0;
  std::size_t totalSamples = 0;
  std::size_t samplesGenerated = 0;
};

struct VolumeComputationResult
{
  std::vector<VolumeComputationRow> rows;
  std::size_t finalSampleCount = 0;
  double volumeEstimate = 0.0;
  std::size_t totalSamplesGenerated = 0;
  std::size_t totalSamplesDeleted = 0;
  double volumeComputationTime = 0.0;
  double samplingTime = 0.0;
  int samplingPrecision = 0;  // GMP decimal digits used for sample bookkeeping
  int latticePrecision = 0;   // b: lattice spacing is 10^-b
};

// Arithmetic backend for the sampling / union-volume bookkeeping.
//   PointRepr (default): the selected walk runs in double, but sampled points
//     are stored and the membership tests (A x <= b) are evaluated in GMP
//     (boost mpf_float) at `samplingPrecision` decimal digits.
//   None (--nogmp):      walk, lattice points and membership tests in double.
//   Full (--fullgmp):    the selected walk runs in GMP at the target precision,
//     then lattice points are stored in GMP.
enum class GmpMode
{
  PointRepr,
  None,
  Full
};

// Random-walk used by volesti to draw (approximately) uniform points from each
// polytope during the union-volume estimation. Selected with --sampler.
//   Billiard (default):    uniform_billiard_walk.hpp -- reflecting trajectories,
//                          mixes fastest for H-polytopes.
//   AcceleratedBilliard:   uniform_accelerated_billiard_walk.hpp.
//   Ball:                  uniform_ball_walk.hpp -- classic ball walk.
//   RDHR / CDHR:           random- / coordinate-directions hit-and-run.
// Under --fullgmp, accelerated billiard falls back to GMP billiard.
enum class SamplerWalk
{
  Billiard,
  AcceleratedBilliard,
  Ball,
  RDHR,
  CDHR
};

// Tunable knobs for the LRA volume engine. A negative walk length means "use
// the volesti default" (volume: 10 + dim/10, sampling: 10).
struct VolumeOptions
{
  int volumeWalkLength = -1;  // --walklen-vol N
  int sampleWalkLength = -1;  // --walklen-samp N
  bool cddSimplify = true;    // cleared by --no-cdd-simp
  // Use the analytic product of side lengths for axis-aligned boxes. Reject
  // any non-box input instead of silently falling back to a volume estimate.
  bool exactBoxVolumes = false;  // --exact-box-volumes
  bool exactParallelotopeVolumes = false;  // --exact-parallelotope-volumes
  // Base seed for the sampling walks and the union-algorithm auxiliary RNG.
  // Set from --seed, but only when the user passes it explicitly, so a run
  // without the flag keeps the historical base seed (123) rather than picking
  // up the CLI default (42).  Each polytope derives its own stream from this
  // base (see lra_volume.cpp), so the walks are not replayed per polytope --
  // which does mean estimates differ from pre-decorrelation runs at any seed.
  unsigned seed = 123;  // --seed N
  SamplerWalk samplerWalk = SamplerWalk::Billiard;  // --sampler
  GmpMode gmpMode = GmpMode::PointRepr;
  // Minimum GMP arithmetic precision, separate from lattice precision b.
  // 0 => auto (b + 32 decimal digits).
  int precision = 0;  // --precision N
  // If non-empty, write each canonicalized polytope to
  // <dumpInePrefix>_cubeN.ine in the cdd/Avis-Fukuda H-representation format
  // (readable by external volume tools such as vinci / lrs). --dump-ine PREFIX.
  std::string dumpInePrefix;
  double epsilon = 0.8;  // --epsilon, 0 < epsilon <= 1
  double delta = 0.2;    // --delta, 0 < delta < 1
};

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const VolumeOptions& options = {},
    const std::function<void(const VolumeComputationRow&)>& onRow = {});

}  // namespace ttc
