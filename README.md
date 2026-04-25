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

Available compsets: `create.newcase --list`. Results appear in `aRawSimuData/` (raw) and `aPlots/` (figures).

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
