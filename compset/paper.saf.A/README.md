# paper.saf.A

Paper-faithful reproduction of Liu et al. (2022) JGR Solid Earth Model A —
Southern San Andreas + San Jacinto fault system.

## Mesh

- `C_mesh = 2` — internal Fortran structured-quad mesh, built by
  `meshgen1.f90` from the per-fault `x{1,2,3}_1.txt` geometry files
  shipped in `comp_input/`.
- No `meshgen.py` step (mesh generation happens inside the binary).

## Loading

- `Rate_direction.txt` (in `comp_input/`) provides the per-fault-node
  loading: column 1 = loading rate (s^-1), column 2 = angle between
  loading direction and fault tangent (deg). This is the smooth
  per-node table recovered from the paper-era local archive — the
  source of the published Model A results.

## Running

```bash
create.newcase --work_dir work/paper.saf.A --compset paper.saf.A
cd work/paper.saf.A
python3 case.setup
bash run.sh
```

`OMP_NUM_THREADS` defaults to 1; override before `bash run.sh` for faster
wall clock. Default `icend = 4000` cycles — Ctrl-C / `kill` whenever you
have enough cycles for the analysis at hand.

## Files

```
user_defined_params.py              # case parameters (icend, friction, viscosity, ...)
comp_input/
  x1_1.txt .. x3_1.txt              # per-fault structured-mesh geometry
  Rate_direction.txt                # per-node loading rate + angle
README.md
```

## Companion compset

`compset/saf.gmsh.lite/` is the C_mesh=3 (gmsh) reproduction of the
same physics on a coarser unstructured mesh, useful for quicker
exploratory runs.
