# EQdyna.2Dcycle

2D finite-element code for physics-based multicycle earthquake dynamics on geometrically complex fault systems.

## Install

```bash
./install.sh             # autodetects OS, installs python deps, builds binary
```

Or explicitly:

```bash
./install.sh -e macos    # macOS  (with python deps)
./install.sh -e ubuntu   # Linux  (with python deps)
./install.sh -m ubuntu   # build only, deps assumed present
```

Sets `EQDYNA2DCYCLEROOT`, adds `bin/` and `scripts/` to `PATH`, and builds `bin/run_eqdyna2d_<VERSION>`.

Dependencies: `gfortran`, Python 3 with `numpy`, `matplotlib`, `xarray`.

## Run

One-step demo (paper.saf.A — Liu et al. 2022 Model A):

```bash
bash example_workflow.sh
```

Manual flow (any compset):

```bash
create.newcase --work_dir work/my_case --compset paper.saf.A
cd work/my_case
python3 case.setup
bash run.sh
```

Available compsets: `create.newcase --list`.

| Compset | Mesh | Fault system |
|---|---|---|
| `paper.saf.A` | C_mesh=2 (Fortran) | Southern San Andreas + San Jacinto — Liu et al. (2022) Model A |
| `saf.gmsh.lite` | C_mesh=3 (gmsh) | the same system on an unstructured mesh, for quicker runs |
| `subei.gmsh.lite` | C_mesh=3 | Subei (atf, dxs, sbt) |
| `gulang.gmsh.lite` | C_mesh=3 | Gulang, ft1..ft5 |
| `xianshuihe.gmsh.lite` | C_mesh=3 | Xianshuihe, 7 faults, per-node loading from GSRM v2.1 |

C_mesh=3 compsets need `python3 meshgen.py` before `bash run.sh`.

## Restart from a previous run

Every cycle the binary refreshes a `binaryop` file holding the full restart state. To extend a finished run from cycle N+1 → M:

```bash
cd work/my_case
# 1. Bump icstart to last_cycle + 1, raise icend.
sed -i 's/par\.icstart, par\.icend = .*/par.icstart, par.icend = N+1, M/' user_defined_params.py
# 2. Regenerate FE_Global.txt + run.sh from the new params.
python3 case.setup
# 3. Launch — binary loads binaryop and resumes at cycle N+1.
bash run.sh
```

Notes:
- For C_mesh=3 (.gmsh.lite) cases `binaryop` already lives in the case dir. For C_mesh=2 (paper.saf.A) it gets moved to `aRawSimuData/` after a finished run — copy it back to the case dir before step 3.
- The new outputs land as `totalop.txt<icstart>` / `cyclelog.txt<icstart>` / `interval.txt<icstart>` (separate from the original `*.txt1`), and run.sh's `rm -f totalop.txt*` only wipes the case dir, so `aRawSimuData/` keeps the previous output safe.
- Post-processing scripts (`plotRuptureDynamics`, `analyze_catalog.py`, `plot_event_slips_overtime_fig4.py`) auto-discover all `totalop.txt*` files via `saf_result_utils.discover_cycle_tags` and stitch them into one continuous catalog.

## Outputs

After a run finishes, `run.sh` moves raw outputs into `aRawSimuData/`:

| File | Content |
|---|---|
| `totalop.txt<N>` | per-cycle stack: shear, normal, slip, slip-rate, rupture-time per fault node, for cycles starting at icstart=N |
| `cyclelog.txt<N>` | one line per cycle: `cycle_id  nucleation_node` |
| `interval.txt<N>` | one value per cycle: interseismic duration (yr) |
| `binaryop` | restart state for resuming with `icstart > 1` |

Plus the input/mesh files used by the run (`vert.txt`, `fac.txt`, `nsmp.txt`, `nsmpGeoPhys.txt`, `meshGeneralInfo.txt`, `Rate_direction.txt`, `FE_*.txt`).

Plots go to `aPlots/`.

## Post-processing

All scripts live in `scripts/` (and are on `PATH` after install). Run from inside a case dir.

| Script | What it does |
|---|---|
| `plotRuptureDynamics` | Per-cycle 4-panel plot (shear, normal, slip, rupture-time vs along-strike). Defaults: only saves cycles with M ≥ `MIN_PLOT_MAGNITUDE` (env var, default 6.5). Set `FORCE_REPLOT=1` to redo existing plots. |
| `CATALOG=1 plotRuptureDynamics` | **Catalog mode**: no figures, ~10× faster. Writes `aPlots/catalog.csv` with columns: `eqId, magnitude, moment_Nm, nuc_x_km, nuc_y_km, nuc_ft, rup_dur_s, peak_slip_m`. |
| `analyze_catalog.py [case_dir] [--mc M] [--mmax M]` | Reads `aPlots/catalog.csv`. Outputs b-value (LSQ on `[Mc, Mmax]` window), magnitude-frequency distribution, magnitude-vs-cycle scatter, nucleation-along-strike-vs-cycle scatter. Saves `aPlots/catalog_analysis.png`. For characteristic-fault MFDs use `--mmax 7.0` to exclude the bump. |
| `plot_event_slips_overtime_fig4.py [case_dir]` | Paper Figure-4 style slip-distribution stacks (slip vs along-strike, vertically offset by event time). `--duration` = window in kyr (default 3); `--threshold` = min event slip (m). |
| `make_paper_figures.py <case_dir>` | One command for the whole paper figure set — catalog, figures 3/4/5/6/9, b-value analysis, per-cycle rupture dynamics. `--only`/`--skip` select stages; each logs to `aPlots/logs/<stage>.log`. |
| `snapshot_case.py <case_dir>` | Snapshots a **running** case for safe post-processing: truncates each segment to whole cycles and merges restart segments (`totalop.txt1` + `totalop.txt73` + ...) into one continuous sequence. Plotting a live or restarted case directly reads half-written cycles and silently sees only one segment. |
| `monitor_runs.sh [-i secs] [-1] [cases...]` | Background poller (default 1200s). Snapshots each named case (or any `work/*/` holding `totalop.txt1`), re-runs the plot suite, and reports events, span, Mmax and the M≥6.5 count. |
| `paleo_site_stats.py` | Table-2 style recurrence and slip statistics at the BF/FM/WW paleoseismic sites (SAF only). |
| `plot_saf_figure6.py`, `plot_saf_figure9.py` | Ports of the published Figure 6 (cumulative moment, magnitude-frequency) and Figure 9 (characteristic-event slip distributions). SAF only. |
| `fetch_published_reference.sh` | Fetches the published Zenodo software and Pangaea results from their DOIs with md5 pinning, instead of vendoring them. |
| `compare_cycle_over_strike.py` | Overlay a chosen cycle's stress/slip/rupture-time curves from multiple cases (e.g. saf.gmsh.lite vs paper.saf.A) for direct comparison. |

Quick example: catalog + analyze on a finished case:

```bash
cd work/my_case
CATALOG=1 plotRuptureDynamics      # writes aPlots/catalog.csv
analyze_catalog.py . --mmax 7.0    # writes aPlots/catalog_analysis.png
```

## Reproducibility

Runs are only reproducible at a **fixed `OMP_NUM_THREADS`**. The OpenMP
reductions in `qdct3` and `hrglss` accumulate under `$OMP ATOMIC`, so the
summation order depends on thread scheduling and results differ in the last
bit between thread counts. Earthquake sequences are chaotic, so that grows
into a visibly different catalogue within about four cycles.

At `OMP_NUM_THREADS=1`, v2.1.0 reproduces the stored SAF references exactly
on both meshing paths — `paper.saf.A` (C_mesh=2) and `saf.gmsh.lite`
(C_mesh=3), max absolute difference 0.0 in `totalop.txt1`. Comparing against
a reference at a different thread count will diverge, and that is not a
regression.

## Mesh checks

C_mesh=3 cases should be checked before a long run:

```bash
python3 scripts/checkMeshQuality.py <case_dir>   # hard-fails (exit 1) on real defects
python3 scripts/plotMeshFaults.py <case_dir>     # per-segment and per-tip views
```

`checkMeshQuality.py` fails on orphaned split nodes, element-count mismatch
against the `.msh`, triangles in a quad mesh, interior angles below 20° or
above 160°, aspect ratio above 10, and degenerate or mixed cells.
`plotMeshFaults.py` renders every fault end to end and every fault **tip** —
the tip view exists because `plotMeshNearFault.py` zooms on midpoints, which
is why every orphaned split node found on this project went unseen.

## Conventions and tests

`PROJECT_RULES.md` records the conventions this project has had to learn the
hard way, most of them from silent failures. `test_system/test_conventions.py`
enforces the mechanical ones:

```bash
python3 test_system/test_conventions.py          # 36 checks
python3 test_system/smoke.py                     # compile + 1-cycle run
python3 -m test_system.test_all                  # full pipeline
```

Releases follow the gates in `PROJECT_RULES.md` (R22-R27): clean tree, docs
touched since the last tag, a non-empty `CHANGELOG.md` body for the version in
`VERSION`, a passing smoke test, then an annotated tag carrying that changelog
section.

## Authors

- **Dunyu Liu** — Institute for Geophysics, Jackson School of Geosciences, The University of Texas at Austin — <dliu@ig.utexas.edu>
- **Benchun Duan** — Center for Tectonophysics, Department of Geology and Geophysics, Texas A&M University — <bduan@tamu.edu>

## Contributors

- **Claude** (Anthropic Claude Code) — development assistance, refactoring, testing, documentation

## Citations

If you use EQdyna.2Dcycle in your research, please cite:

- Duan, B., & Oglesby, D. D. (2006). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults. *Journal of Geophysical Research*, 111(B5), B05309. <https://doi.org/10.1029/2005JB004138>
- Liu, D., Duan, B., Scharer, K., & Yule, D. (2022). Observation-constrained multicycle dynamic models of the southern San Andreas and the northern San Jacinto faults: Addressing complexity in paleoearthquake extent and recurrence with realistic 2D fault geometry. *Journal of Geophysical Research: Solid Earth*, 127(2), e2021JB023420. <https://doi.org/10.1029/2021JB023420>

## License

MIT — see [LICENSE](LICENSE).
