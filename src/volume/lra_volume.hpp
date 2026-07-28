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
};

// Arithmetic backend for the sampling / union-volume bookkeeping.
//   PointRepr (default): the billiard walk runs in double, but sampled points
//     are stored and the membership tests (A x <= b) are evaluated in GMP
//     (boost mpf_float) at `samplingPrecision` decimal digits.
//   None (--nogmp):      the original all-double path (no GMP, no truncation).
//   Full (--fullgmp):    the billiard walk itself runs in GMP at the target
//     precision (patched volesti), then points are stored in GMP.
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
// Note: --fullgmp only implements the billiard walk; other choices fall back to
// the billiard walk (with a warning) under --fullgmp.
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
  SamplerWalk samplerWalk = SamplerWalk::Billiard;  // --sampler
  GmpMode gmpMode = GmpMode::PointRepr;
  // GMP precision (decimal digits) for sampling. 0 => auto, derived from the
  // polytopes via get_precision_from_cubes (src/cube_processor_nondis.py).
  int precision = 0;  // --precision N
  // If non-empty, write each canonicalized polytope to
  // <dumpInePrefix>_cubeN.ine in the cdd/Avis-Fukuda H-representation format
  // (readable by external volume tools such as vinci / lrs). --dump-ine PREFIX.
  std::string dumpInePrefix;
};

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const VolumeOptions& options = {},
    const std::function<void(const VolumeComputationRow&)>& onRow = {});

}  // namespace ttc

