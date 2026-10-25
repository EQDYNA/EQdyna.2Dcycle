# EQdyna.2Dcycle

2D finite-element code for physics-based multicycle earthquake dynamics on
geometrically complex fault systems.

**Full user guide:** <https://eqdyna.github.io/EQdyna.2Dcycle/> — install
details, every compset, parameters, outputs, post-processing, source time
functions, testing/reproducibility, and troubleshooting.

## Install

```bash
git clone https://github.com/EQDYNA/EQdyna.2Dcycle.git
cd EQdyna.2Dcycle
./install.sh             # autodetects OS, installs python deps, builds the binary
```

Or explicitly: `./install.sh -e macos`, `./install.sh -e ubuntu` (build +
python deps), or `./install.sh -m <machine>` (build only, deps assumed
present). This sets `EQDYNA2DCYCLEROOT`, adds `bin/` and `scripts/` to
`PATH`, and builds `bin/run_eqdyna2d_<VERSION>`.

Dependencies: `gfortran` (with OpenMP), Python 3 with `numpy`, `matplotlib`,
`xarray`. See [Getting started](https://eqdyna.github.io/EQdyna.2Dcycle/getting-started/)
for gmsh/meshio (needed only for the `*.gmsh.lite` compsets) and other
optional packages.

## Run one example

```bash
bash example_workflow.sh
```

This builds the code, creates a case from the `paper.saf.A` compset
(reproducing Liu et al. 2022 JGR Solid Earth Model A — southern San Andreas
+ northern San Jacinto), and launches it in the background. It runs 4000
earthquake cycles by default; `Ctrl-C` once you have enough cycles for what
you're testing, or lower `par.icend` in
`work/paper.saf.A.demo/user_defined_params.py` first.

Equivalently, step by step:

```bash
create.newcase --work_dir work/my_case --compset paper.saf.A
cd work/my_case
python3 case.setup   # writes FE_*.txt + run.sh from user_defined_params.py
bash run.sh           # runs in the background; logs to run_<timestamp>.log
```

`create.newcase --list` shows all compsets. The `*.gmsh.lite` compsets need
`python3 meshgen.py` before `bash run.sh`; see
[Compsets](https://eqdyna.github.io/EQdyna.2Dcycle/compsets/).

## Outputs

When a run finishes, `run.sh` moves the raw results — `totalop.txt<N>`,
`cyclelog.txt<N>`, `interval.txt<N>`, `binaryop` — into `aRawSimuData/` in
the case directory, and plots go to `aPlots/`. See
[Output files](https://eqdyna.github.io/EQdyna.2Dcycle/outputs/) for the
full format and [Post-processing](https://eqdyna.github.io/EQdyna.2Dcycle/post-processing/)
for the plotting and analysis scripts (`plotRuptureDynamics`,
`analyze_catalog.py`, `make_paper_figures.py`, ...).

## Getting help

See [Troubleshooting](https://eqdyna.github.io/EQdyna.2Dcycle/troubleshooting/),
or open a [GitHub issue](https://github.com/EQDYNA/EQdyna.2Dcycle/issues).
For development conventions and release process, see `PROJECT_RULES.md` and
`CHANGELOG.md`.

## Authors

- **Dunyu Liu** — Institute for Geophysics, Jackson School of Geosciences, The University of Texas at Austin — <dliu@ig.utexas.edu>
- **Benchun Duan** — Center for Tectonophysics, Department of Geology and Geophysics, Texas A&M University — <bduan@tamu.edu>

## Citing

If you use EQdyna.2Dcycle in your research, please cite:

- Duan, B., & Oglesby, D. D. (2006). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults. *Journal of Geophysical Research*, 111(B5), B05309. <https://doi.org/10.1029/2005JB004138>
- Liu, D., Duan, B., Scharer, K., & Yule, D. (2022). Observation-constrained multicycle dynamic models of the southern San Andreas and the northern San Jacinto faults: Addressing complexity in paleoearthquake extent and recurrence with realistic 2D fault geometry. *Journal of Geophysical Research: Solid Earth*, 127(2), e2021JB023420. <https://doi.org/10.1029/2021JB023420>

## License

MIT — see [LICENSE](LICENSE).
