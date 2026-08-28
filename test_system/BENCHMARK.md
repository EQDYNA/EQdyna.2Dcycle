# Reference benchmarks — EQdyna.2Dcycle

Wall-clock timings + summary catalogs for `icend = 5000` runs of the three
`*.gmsh.lite` compsets, captured as a regression reference. Run with
`EQdyna.2Dcycle v2.0.7-rc7` (binary `run_eqdyna2d_2.0.7-rc7`).

## Hardware

| | |
|---|---|
| Host | `cotopaxi` |
| OS | Ubuntu 22.04.5 LTS, kernel Linux 6.8.0-100-generic x86_64 |
| CPU | 2× AMD EPYC 7532 (32 cores/socket, 1 thread/core — 64 cores total, no SMT). 1.5–2.4 GHz. |
| Cache | L1d 64×32 KiB, L1i 64×32 KiB, L2 64×512 KiB, L3 32×16 MiB (512 MiB total) |
| NUMA | 8 nodes × 8 cores |
| Memory | 1 TiB |

## Toolchain / build

| | |
|---|---|
| Compiler | GNU Fortran 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3) |
| Python | 3.12.4 |
| Build flags | `-O3 -march=native -funroll-loops -ffree-line-length-none -std=legacy -fopenmp -cpp -DEQDYNA_VERSION=...` |
| Run env | `OMP_NUM_THREADS=1` (default) |

All runs were launched in parallel via `bash run.sh` with the default
serial-Fortran path (one core per case). The 64-core box is otherwise
idle, so per-cycle timing is not affected by competing load.

## Per-compset wall time

Each compset was first run cold from `icstart=1` to its shipped `icend`,
then restarted from `binaryop` to `icend=5000` (see README "Restart from
a previous run"). The two segments are reported separately because their
event distributions differ.

| Compset | Segment | Cycles | Wall time | min / 1000 cycles |
|---|---|---:|---:|---:|
| saf.gmsh.lite     | 1..3000   | 3000 | 129 min  | 43.0 |
| saf.gmsh.lite     | 3001..5000 | 2000 | 72 min   | 36.0 |
| **saf.gmsh.lite (total)** | **1..5000** | **5000** | **201 min** | **40.2** |
| subei.gmsh.lite   | 1..1000   | 1000 | 52 min   | 52.0 |
| subei.gmsh.lite   | 1001..5000 | 4000 | 424 min  | 106.0 |
| **subei.gmsh.lite (total)** | **1..5000** | **5000** | **476 min** | **95.2** |
| gulang.gmsh.lite  | 1..4000   | 4000 | 146 min  | 36.5 |
| gulang.gmsh.lite  | 4001..5000 | 1000 | 32 min   | 32.0 |
| **gulang.gmsh.lite (total)** | **1..5000** | **5000** | **178 min** | **35.6** |

Notes on per-1000-cycle variability:
- The interseismic-loading step is fast; almost all wall time is in the
  dynamic-rupture solve. So per-1000-cycle timing is driven by how many
  large (M ≳ 6) events fall in that range.
- subei's later cycles (1001..5000) are 2× slower per-1000 than its
  early cycles (1..1000) because the recurrent M ~ 7 events on the long
  ATF cluster densely after warm-up.
- gulang and saf are roughly stationary across the cold + restart
  segments (~36 vs 32 min/1000, ~43 vs 36 min/1000 respectively).

## Catalog statistics (cold segment, icstart=1)

Captured from `aPlots/catalog.csv` written by
`CATALOG=1 plotRuptureDynamics`, with b-value from
`analyze_catalog.py . --mmax 7.0` (LSQ fit on `[Mc, 7.0]`).

| Compset | Cycle range | events | M range | b-value (Mc..7.0) |
|---|---|---:|---|---:|
| saf.gmsh.lite    | 1..3000 | 3000 | 4.87 .. 8.02 | 0.86 ± 0.07 |
| subei.gmsh.lite  | 1..1000 | 1000 | 4.77 .. 7.34 | 0.89 ± 0.07 |
| gulang.gmsh.lite | 1..4000 | 4000 | 2.03 .. 7.32 | 0.66 ± 0.03 |

(`plotRuptureDynamics` catalog mode reads `totalop.txt<icstart>` from
`FE_Global.txt`, so it captures one restart segment at a time. Merged
cross-segment statistics for the full 1..5000 run are deferred — the
underlying multi-segment loader in `saf_result_utils` already handles
the stitching for `analyze_catalog.py` and `plot_event_slips_overtime_fig4.py`.)

## paper.saf.A reference

The `paper.saf.A` run (C_mesh=2, dxy=200 m — finer mesh than the lite
compsets) completed `icstart=1 .. icend=4000` in `work/paper.saf.A.demo`.

| | |
|---|---|
| Binary | `run_eqdyna2d_2.0.7-rc6` |
| Started | 2026-04-25 08:41:19 CDT |
| Last cycle written | 2026-05-10 07:38:52 CDT |
| Cycles | 4000 |
| Wall time | 14 d 22 h 58 m (21 478 min) |
| min / 1000 cycles | 5369 |

That is ~134× the per-cycle cost of `saf.gmsh.lite` (40.2 min/1000), consistent
with the 2× finer mesh (200 m vs 400 m) plus the C_mesh=2 structured-mesh
domain being larger: per-cycle wall time is dominated by `qdct3` in the
dynamic-rupture solve (~2 s per simulated millisecond vs ~0.3 s for the lite
compsets).

Catalog (`CATALOG=1 plotRuptureDynamics`, b-value from
`analyze_catalog.py . --mmax 7.0`):

| Compset | Cycle range | events | M range | Mc | b-value (Mc..7.0) |
|---|---|---:|---|---:|---:|
| paper.saf.A | 1..4000 | 4000 | 3.93 .. 8.23 | 3.90 | 0.585 ± 0.055 |

The lower b-value than `saf.gmsh.lite` (0.86) reflects the finer mesh
resolving more of the small-event tail *and* the characteristic M ~ 8 bump
sitting just above the LSQ window — the fit is on `[3.9, 7.0]` in both cases.

Note: the run's nohup wrapper was lost before the binary exited, so `run.sh`
never executed its post-run steps. They were completed by hand afterwards
(`plotRuptureDynamics`, then `mv totalop.txt* cyclelog.txt* interval.txt*
binaryop aRawSimuData/`), so the case is now in the normal post-run layout and
`analyze_catalog.py` / `plot_event_slips_overtime_fig4.py` resolve it via the
`aRawSimuData/` search path. The cycle data is complete — `cyclelog.txt1`
reaches 4000. Wall time above is `Job started` in the log to the mtime of the
last-written output, since the log carries no `Job finished` timestamp.

## How to reproduce

```bash
source install.sh
for c in saf.gmsh.lite subei.gmsh.lite gulang.gmsh.lite; do
  d=work/$c.bench
  create.newcase --work_dir $d --compset $c --force
  ( cd $d && python3 case.setup && bash run.sh )
done
```

After the cold runs finish (icend per compset), restart each to 5000 via
the README "Restart from a previous run" recipe. Wall times above were
captured from `run_*.log` `Job started` / `Job finished` markers.
