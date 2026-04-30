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
| `monitor_runs.sh [seconds]` | Background poller (default 600s). Tails progress for all `paper.saf.A.*` and `saflite` cases under `work/`, re-runs Figure 4 + rupture-dynamics plots each tick. |
| `compare_cycle_over_strike.py` | Overlay a chosen cycle's stress/slip/rupture-time curves from multiple cases (e.g. saf.gmsh.lite vs paper.saf.A) for direct comparison. |

Quick example: catalog + analyze on a finished case:

```bash
cd work/my_case
CATALOG=1 plotRuptureDynamics      # writes aPlots/catalog.csv
analyze_catalog.py . --mmax 7.0    # writes aPlots/catalog_analysis.png
```

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
