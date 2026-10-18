# Combined volume and sampling error benchmarks

These `VOLANDSAMP_*` inputs are two-sided axis-spoke boxes with exactly known
true volumes, but TTC runs its ordinary stochastic body-volume estimator.
The individual body estimate can be wrong, so a union error here does not
isolate a sampling failure. The [`manifest.csv`](manifest.csv),
[`results.csv`](results.csv), and [`benchmarks/`](benchmarks/) contain the
geometries, exact comparisons, and TTC logs.

There are six benchmark geometries and 20 saved TTC runs. The case
`VOLANDSAMP_m8_w20_v1` fails at sampling walk **500**: TTC returns 408.582
against exact union 1168, a 65.02% error. Its implied overlap at walk 500
is 8.86%, close to the exact 8.75%; its body-volume estimates are the main
source of that persistent error. The `v50` cases compare 8, 10, 12, 16, and
20 bodies at one fixed width and volume walk. Their figures share axes across
dimensions and label only body count. All figures are under [`figures/`](figures/)
as PDF, PGF, and PNG.

Reproduce from the repository root:

```bash
python3 example/kakeya/vol_and_samp/generate.py
python3 example/kakeya/vol_and_samp/run.py --missing
python3 example/kakeya/vol_and_samp/analyze.py
python3 example/kakeya/vol_and_samp/plot.py
```

## 20 target walks × 5 ordinary-volume cases

[`target_manifest.csv`](target_manifest.csv) contains 100 two-body cases at
the same 20 target walks and five width variants as `../samp_only/`.
Through walk 127, the input geometry also matches exactly. Above 127,
the ordinary-volume cases use thin denominators 1009–1019 instead of
`20 × target walk`: the original, thinner version made the volume oracle
error or stall before sampling. Those original inputs and logs are archived
at `example/kakeya_previous/target_v1_vol_and_samp_high/`.
TTC uses its ordinary volesti volume estimate with `--walklen-vol 1` and
the ball sampler at the target walk length. The manifest records analytic
true volumes, and
[`target_results.csv`](target_results.csv) records measured body-volume and
union errors only when TTC finishes. The
[`target_union_error_measured.pdf`](figures/target_union_error_measured.pdf),
[`target_run_coverage.pdf`](figures/target_run_coverage.pdf), and
[`target_overlap_mismatch_measured.pdf`](figures/target_overlap_mismatch_measured.pdf)
figures are also saved as PGF and PNG, as is
[`target_projection_x0_x1.pdf`](figures/target_projection_x0_x1.pdf).
The overlap plot compares sampled deletion with the analytic body
intersection. Unrun targets are kept out of the error graph.
With `--walklen-vol 1`, the archived thinner targets at walks 234 and 428
returned invalid body volumes before sampling, and one target-785 probe
timed out in volume computation after 120 seconds. They are not counted as
measured union errors. TTC prints body estimates with six fixed decimals,
so very small positive estimates may appear as `0.000000`; body error and
apparent overlap are left blank when those printed values round to zero.

```bash
python3 example/kakeya/vol_and_samp/run_targets.py --targets 127 --missing
python3 example/kakeya/vol_and_samp/analyze_targets.py
python3 example/kakeya/vol_and_samp/plot_targets.py
```
