# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EQdyna.2Dcycle is a 2D finite element method (FEM) code for simulating physics-based multicycle earthquake dynamics on geometrically complex fault systems. The codebase consists of:

- **Core simulation engine**: Fortran 90 source code in `src/` 
- **Case management**: Python scripts for creating and configuring simulation cases
- **Post-processing**: Python utilities for plotting and analysis
- **Test suites**: Reference results and validation tools

## Build and Installation

### Install dependencies and build executable:
```bash
# For macOS (sets up environment and builds)
./install.sh -e macos

# For Ubuntu/Linux
./install.sh -m ubuntu

# Just build (assumes dependencies exist)
./install.sh -m macos
```

### Manual build:
```bash
cd src/
make           # produces run_eqdyna2d_<VERSION> using the top-level VERSION file
cd ..
mkdir -p bin
cp src/run_eqdyna2d_* bin/
```

The root `VERSION` file (single-line, e.g. `2.0.3`) is read by `src/makefile`;
the binary is named `run_eqdyna2d_$(VERSION)`. Bump the version there and
rebuild to produce a new binary alongside old ones.

### Environment setup:
```bash
# Source to set environment variables
source install.sh
```

This sets `EQDYNA2DCYCLEROOT` and adds `bin/` and `scripts/` to PATH.

## Running Simulations

### Create a new case:
```bash
create.newcase --work_dir work/my_case --compset paper.saf.A          # named flags (preferred)
create.newcase work/my_case paper.saf.A                                # legacy positional still works
create.newcase --list                                                  # list available compsets
create.newcase --work_dir work/my_case --compset paper.saf.A --force   # overwrite existing
```

`create.newcase` fails by default if the destination exists; pass `--force`
to wipe and recreate. It copies `compset/<name>/*`, `scripts/*`, and `bin/*`
into the case dir (skipping generated files like FE_*.txt and run.sh).

### Configure and run:
```bash
cd work/my_case
# Edit user_defined_params.py if needed
python3 case.setup   # writes FE_*.txt + run.sh from user_defined_params.py
# For C_mesh=3 (gmsh) cases:
python3 meshgen.py   # produces fem_mesh_output/ and nsmpGeoPhys.txt
bash run.sh          # launches the simulation in background via nohup
```

The generated `run.sh` exports `GFORTRAN_UNBUFFERED_ALL=1` (flushes Fortran
stdout immediately for real-time log tailing) and
`OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}` (default serial; can be overridden in the caller's
env). Output goes to a timestamped `run_YYYYMMDD_HHMMSS.log`; raw results
(`totalop.txt*`, `cyclelog.txt*`, `interval.txt*`, `binaryop`) are moved to
`aRawSimuData/` after the binary exits.

## Testing

### Run automated test suite:
```bash
python3 -m test_system.test_all  # Complete automated testing pipeline
```

This runs:
- Fresh compilation from source
- Case creation with organized folder structure  
- Mesh generation and setup
- Full simulation execution
- Post-processing and visualization
- Verification against reference results

### Manual verification:
```bash
python3 test_system/verify.test.py  # Compare results against reference data
```

## Key File Types and Architecture

### Core Fortran Source (`src/`)
- `eqdyna2d.f90`: Main program entry point
- `globalvar.f90`: Global variables and data structures
- `driver.f90`: Main simulation driver
- `faulting.f90`, `fric.f90`: Fault mechanics and friction laws
- `meshgen.f90`, `meshgen1.f90`: Mesh generation
- `qdct*.f90`: Finite element routines (2D/3D quadrilateral elements)
- `contm*.f90`: Contact mechanics
- `interstress.f90`: Interseismic stress calculations

### Case Configuration
- `compset/paper.saf.A/`: paper-faithful Liu et al. (2022) Model A reproduction (C_mesh=2). Uses the smooth per-node Rate_direction.txt recovered from the paper-era local archive.
- `compset/saf.gmsh.lite/`: gmsh-meshed (C_mesh=3) paper-A reproduction at coarser dxy=400m. Self-contained — ships with paper.saf.A's Rate_direction.txt + x{1,2,3}_1.txt; meshgen.py interpolates per-node loading onto gmsh fault nodes via natural-cubic-spline analytical tangent. First cycle: 473 yr / M5.7 (paper-A: 471 yr / M5.74).
- `compset/subei.gmsh.lite/`: gmsh-meshed (C_mesh=3) Subei fault system (atf, dxs, sbt). Same framework as saf.gmsh.lite but uses the default uniform-x loading mode (`uniformXLoadingAngle`) since no paper Rate_direction.txt is available. Multi-surface decomposition handles the connected fault junctions.
- `compset/xianshuihe.gmsh.lite/`: Xianshuihe fault (鲜水河断裂), eastern Tibet. **Geometry input only — not yet runnable.** Ships the 1:1M Chinese active-fault KML, `plot_fault_trace.py` (parses the 9 polylines into **8 faults**, ft1..ft8 SE→NW via `MERGE_GROUPS` — seg1+seg2 merged since they meet at 0.26 km and seg1 alone is only 15 km; all faults vertical, `FT_DIP_DEG=90`, `FT_TYPE=1`) and `map_step_overs.py` (junction and splay-strand step-overs, classified releasing/restraining with the sinistral convention). Loading plan from Qiao et al. (2022, EPSL 596, 117799) is in its README.
- `compset/gulang.gmsh.lite/`: gmsh-meshed (C_mesh=3) Gulang fault system (ft1..ft5, all left-strike vertical). Uses uniform-x loading. Ships with `preprocess_control_points.py` that resamples the dense raw `.gmt` digitisation (412 control pts, median ~30 m spacing) down to a SAF-like 39-point set at ≥3 km arc-length spacing via per-fault smoothing splines; originals preserved as `*.gmt.txt.orig`. ftVis tuned to 5e21 Pa·s (1e21 too low for nucleation, 8.4e21 gives tensile normal at large strikes).
- `user_defined_params.py`: Simulation parameters (inherits from defaultParameters.py)
- `userDefinedFaultSysGeoPhys.py`: Fault geometry and physics
- `FE_*.txt`: Input files generated by case.setup
- For the C_mesh=2 `Rate_direction.txt` file format, units, and how it
  feeds into paper eqs (1)–(3), see `archive/RATE_DIRECTION_TXT.md` (untracked local reference).

### Python Utilities (`scripts/`)
- `case.setup`: Generates FE_*.txt input files + run.sh from user_defined_params.py
- `create.newcase`: Creates new case from a compset (supports `--work_dir`, `--compset`, `--force`, `--list` flags; legacy positional form still accepted)
- `plotRuptureDynamics`: Per-cycle shear/normal/slip/rupture-time plots from totalop.txt
- `plot_event_slips_overtime_fig4.py`: Figure-4-style slip-distribution stacks. Default tstart auto-populates 0, D, 2D, … covering total simulated kyr (window duration D = 3 kyr by default).
- `plot_saf_figure3.py`: Figure-3-style long-term slip-rate comparisons.
- `compare_cycle_over_strike.py`: Overlay a chosen cycle's stress/slip/rupture-time curves from multiple cases (plus Pangaea reference) for direct comparison.
- `plotRuptureDynamics`: per-cycle 4-panel plot. `MIN_PLOT_MAGNITUDE` env var (default 6.5) gates which cycles are saved as PNGs; `FORCE_REPLOT=1` re-renders existing plots; **`CATALOG=1`** runs in catalog-only mode (no figures, ~10× faster) and writes `aPlots/catalog.csv` with eqId, magnitude, moment, nucleation x/y, rupture duration, peak slip per cycle.
- `analyze_catalog.py`: reads `aPlots/catalog.csv`, computes b-value (LSQ on user-windowed magnitude range; characteristic-fault MFDs need `--mmax 7.0` to exclude the bump), MFD, magnitude-vs-cycle, nucleation-vs-cycle. Saves `aPlots/catalog_analysis.png`.
- `monitor_runs.sh`: periodic status + Figure 4 + rupture-dynamics re-plot for active `paper.saf.A.*` and `saflite` cases. Defaults to 600s polling; pass seconds as first arg to override.
- `plot_saf_figure6.py`: Figure-6-style cumulative moment release, magnitude histogram, and magnitudes >6.6 against the Scharer & Yule (2020) prehistoric/historical datasets.
- `plot_saf_figure9.py`: Figure-9-style characteristic-event slip distributions. Selects events by **rupture footprint**, not by cycle id — the published MATLAB hard-codes ids from its own run, which mean nothing in another sequence. `--published` restores the paper's exact panels; `--cycles` takes an explicit list.
- `paleo_site_stats.py`: Table-2-style recurrence-interval and slip statistics at the BF/FM/WW paleoseismic sites. Parity-gated on running the published MATLAB against the published Pangaea output.
- `make_paper_figures.py`: **one command for the full paper figure set on any case** — `python3 make_paper_figures.py <case_dir>`. Stages: catalog, figure3/4/5/6/9, analysis, rupture. `--skip rupture` for the cheap summary set, `--only <stages>`, `--list`. Each stage logs to `aPlots/logs/<stage>.log`; the expensive per-cycle stage runs last.
- `fetch_published_reference.sh`: fetches the published Zenodo software (10.5281/zenodo.5823021) and Pangaea results (10.1594/PANGAEA.940262) from their DOIs with md5 pinning, instead of vendoring 615 MB of immutable data.
- `saf_result_utils.py`: Loader helpers used by the plot_saf_* scripts.
- `meshGenLib.py`: Mesh generation utilities (C_mesh=3)
- `defaultParameters.py`: Default simulation parameters

### Output Files
- `totalop.txt*`: Main simulation output files (stress, slip, rupture time)
- `cyclelog.txt*`: Earthquake cycle logging
- `interval.txt*`: Inter-event intervals
- `aRawSimuData/`: Directory containing raw simulation results
- `aPlots/`: Generated plots and figures

## C_mesh Branching Architecture

The `C_mesh` parameter in `user_defined_params.py` controls the mesh and loading workflow:

### C_mesh=2 — Internal Fortran meshing (structured quad mesh)
- Mesh generated internally by `meshgen1.f90` from `x*_1.txt` geometry files
- Loading specified via **`Rate_direction.txt`**: `maxftnode × ntotft` entries, 2 columns per row:
  - Col 1: loading rate (s⁻¹)
  - Col 2: angle between loading direction and fault tangent (degrees)
- Static files (`x*_1.txt`, `Rate_direction.txt`) live in `compset/<case>/comp_input/` and are copied by `case.setup`
- Reference case: `compset/paper.saf.A/` (reproduces Liu et al. 2022 Model A)

### C_mesh=3 — External Python/gmsh meshing (unstructured quad mesh)
- Mesh generated by `meshgen.py` using gmsh; outputs go to `fem_mesh_output/`
- Loading specified via **`nsmpGeoPhys.txt`**: `maxftnode × ntotft` entries, 9 columns per row:
  1. `tx` — fault unit tangent x-component
  2. `ty` — fault unit tangent y-component
  3. `len` — segment length (m)
  4. `ftType` — fault type (1: left-strike, -1: right-strike, 2: thrust, -2: normal)
  5. `ftDip` — fault dip (degrees; 90 = strike-slip)
  6. `ftLoadMaxShear` — far-field shear strain rate (s⁻¹; typically 1.427e-14)
  7. `ftLoadAngle` — angle between loading and fault tangent (degrees; -999 = auto)
  8. `ftLoadWt` — loading weight (default = 450; controls approach timescale only)
  9. `ftVis` — viscosity (Pa·s; typically 6e21)
- `nsmpGeoPhys.txt` is written by `meshgen.py` via `meshGenLib.py`
- Reference cases: `compset/saf.gmsh.lite/` (paper-faithful), `compset/subei.gmsh.lite/`, `compset/gulang.gmsh.lite/`

### ftLoadWt and ftLoadAngle conventions (C_mesh=3)
- **ftLoadWt=450** is the default (full loading). Formula from `strbld.f90`:
  `ant = ftVis * 450 / ftLoadWt`, `rs = ftLoadMaxShear / 450 * ftLoadWt * cos(2θ) * ant`
  → ftLoadWt cancels in steady-state; it only controls the approach timescale.
  Fractional loading examples: `ftLoadWt = 0.8 * 450` (80%), `0.1 * 450` (10%).
- **ftLoadAngle=-999** (auto): angle computed as `θ = atan(ty/tx)` — fault tangent from x-axis,
  assuming tectonic loading is along the general strike direction (x-axis).
- θ is clamped to 45° maximum before computing shear/normal loading components.

### Source directories
- `src/` — the single development baseline (v2.0.3 code path). Has C_mesh=2 (paper/Fortran mesh) and C_mesh=3 (gmsh) branches.
- `archive/src.v2.1.0/` — archived v2.1.0 (gmsh-only branch)
- `archive/published/source_local_recovered/code/` — preserved paper-era v2.0.2 snapshot; read-only local reference (archive/ is untracked).

### Interseismic viscosity (C_mesh=2, `interstress.f90`)
Per-node viscosity varies inversely with local loading rate (Liu et al. 2022 / Duan & Oglesby 2005):
```
ant_i = ant0 * str / rd(1,i)    (rd in rad/s, str in rad/s from FE_Global line 23)
rs_i  =  rd(1,i) * cos(2θ) * ant_i     → simplifies to str * cos(2θ) * ant0
rn_i  = -rd(1,i) * sin(2θ) * ant_i
```
At the fault node with peak `rd(1)`, `ant` hits its minimum `η_min = ant0 * str / rd_max`, which the paper reports in Table 1 (Model A: 7.0×10²¹ Pa·s with ant0=8.4e21 and rd_max≈540 nrad/yr). Fallback `ant = ant0` for nodes with `rd(1)=0` (fault endpoints).

### Performance notes (driver hot loop)
- `alhs_inv = 1/alhs` is precomputed once after `qdct2` (in `eqdyna2d.f90`). The driver uses whole-array `brhs = brhs * alhs_inv` — vectorizes to AVX2 `vmulpd`, avoiding per-timestep divisions.
- `qdct3.f90` skips the `zeroal` check since `rdampm(1)=0` is enforced in `eqdyna2d.f90`; `formma` is always false.
- OpenMP parallelized: `qdct3` (with `$OMP ATOMIC` on brhs assembly), `hrglss` (same pattern), the velocity/displacement update in `driver.f90`, and the `brhs = brhs*alhs_inv` multiply. Thread count is controlled by `OMP_NUM_THREADS` in run.sh. Scaling on 4 threads: ~2.8–3× (qdct3 limited by atomic contention; hrglss near-linear).
- Per-op timing (`t_vd`, `t_qdct3`, `t_hrglss`, `t_faulting`, `t_brhs`) is printed every `print_every` timesteps (default 300 = every 3 simulated seconds) in the run log, alongside a date/time stamp.

### Magnitude convention (important when comparing to the paper)
`plotRuptureDynamics` computes moment with the model's own constants,
`mu = rou*vs^2 = 2670*3464^2 = 3.204e10 Pa` at a 22 km seismogenic depth.
The paper's Figures 6 and 9 instead use `mu = 3500^2*3000 = 3.675e10 Pa`
at the same depth, so **every magnitude and b-value from `catalog.csv`
sits 0.04 Mw below the paper's published scale**. `plot_saf_figure6.py`
and `plot_saf_figure9.py` use the paper's constants and are therefore NOT
directly comparable to `catalog.csv`. The paper's constant is itself
inconsistent with its own model parameters (rou=2670, vs=3464); making
this selectable is planned.

### Output files
- `meshGeneralInfo.txt` — mesh summary (always written; replaces `Mesh_general_info.txt`)
- `nsmpTanLen.txt` — same content as `nsmpnv.txt` (fault tangent + segment lengths)
- `totalop.txt{icstart}` — all cycles stacked (ncycles × totftnode rows, 5 cols per node)
- Fortran unit 6 = stdout; output file uses **unit 61** to avoid stdout conflict

## Code Architecture

### Simulation Workflow
1. **Case Setup**: Configure parameters in `user_defined_params.py`
2. **Input Generation**: Run `case.setup` to create FE_*.txt files
3. **Mesh Generation**: Fortran code (C_mesh=2) or `meshgen.py` (C_mesh=3)
4. **Fault Loading**: Apply tectonic loading and build up stress
5. **Dynamic Rupture**: Simulate earthquake rupture dynamics
6. **Post-Processing**: Generate plots and analyze results

### Key Components
- **Friction Laws**: Multiple friction models (friclaw parameter)
- **Mesh Types**: Structured (C_mesh=2) or unstructured (C_mesh=3) quadrilateral elements
- **Fault Geometry**: Complex multi-segment fault systems
- **Loading**: Quasi-static interseismic loading with viscosity
- **Output**: Time series of stress, slip, and rupture properties

### Parameter Management
Parameters flow from `defaultParameters.py` → `user_defined_params.py` → `case.setup` → FE_*.txt input files. The `par` object contains all simulation parameters including mesh resolution, friction coefficients, material properties, and fault geometry.

## Dependencies
- **Fortran**: gfortran compiler
- **Python 3**: numpy, matplotlib, xarray
- **Optional**: gmsh, meshio, nbconvert (for advanced meshing)

## Translation Notes
When translating MATLAB code to Python, maintain identical algorithms, file names, and code structure. Focus on line-by-line equivalence while adapting syntax for Python conventions.