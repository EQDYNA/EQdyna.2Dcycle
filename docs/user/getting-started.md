# Getting started

## Requirements

| | |
|---|---|
| Fortran | `gfortran` with OpenMP |
| Python 3 | `numpy`, `matplotlib`, `xarray` |
| gmsh meshing (C_mesh = 3 cases only) | `gmsh` Python package, `meshio` |
| Optional | `scipy` (some analysis scripts), `netCDF4` (NetCDF export) |

## Install and build

```bash
git clone https://github.com/EQDYNA/EQdyna.2Dcycle.git
cd EQdyna.2Dcycle
./install.sh             # detects the OS, installs Python dependencies, builds
```

Or choose explicitly:

```bash
./install.sh -e ubuntu   # Linux, with Python dependencies
./install.sh -e macos    # macOS, with Python dependencies
./install.sh -m ubuntu   # build only; dependencies already present
```

The build produces `bin/run_eqdyna2d_<VERSION>`, where `<VERSION>` is the
contents of the top-level `VERSION` file. To use the tools in a new shell,
source the installer once:

```bash
source install.sh        # sets EQDYNA2DCYCLEROOT, puts bin/ and scripts/ on PATH
```

To rebuild by hand after changing the Fortran:

```bash
cd src && make && cd ..
mkdir -p bin && cp src/run_eqdyna2d_* bin/
```

## A first run

`paper.saf.A` is Model A of Liu et al. (2022): the southern San Andreas and
northern San Jacinto faults. A short run:

```bash
create.newcase --work_dir work/first --compset paper.saf.A
cd work/first
sed -i 's/^par.icstart, par.icend = .*/par.icstart, par.icend = 1, 3/' user_defined_params.py
python3 case.setup
bash run.sh
```

`run.sh` runs in the background and writes a timestamped log,
`run_YYYYMMDD_HHMMSS.log`. When it finishes, the raw results are in
`aRawSimuData/`. The first three interseismic intervals should be 471, 33
and 19 years:

```bash
cat aRawSimuData/interval.txt1
```

For the whole Model A demonstration in one step, `bash example_workflow.sh`
from the repository root.
