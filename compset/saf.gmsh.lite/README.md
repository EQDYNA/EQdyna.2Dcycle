# saf.gmsh.lite

Coarser-resolution gmsh-meshed (C_mesh=3) reproduction of the Southern San
Andreas Fault system from Liu et al. (2022) Model A. Companion to
`paper.saf.A` (C_mesh=2 reference). Standard build/run flow: repo root
`README.md` § Run (`create.newcase` → `python3 meshgen.py` → `case.setup` →
`bash run.sh`).

Three independent right-lateral strike-slip faults (ftType=-1, ftDip=90):
`ssaf` — Southern San Andreas Fault (80 control points, single trace); `sjfn`
— San Jacinto Fault North Strand (15 control points); `sjfs` — San Jacinto
Fault South Strand (10 control points). Control points:
`user_fault_geometry_input/{ssaf,sjfn,sjfs}.gmt.txt` (cols 1–2 = x, y in km;
cols 3–4 unused, superseded by paper interpolation).

`meshgen.py` builds a natural cubic spline y(x) per fault (mirrors
`meshgen1.f90`), meshes with gmsh (`dx=0.4 km`), and for every gmsh fault node
interpolates per-node loading from `paper.saf.A/Rate_direction.txt` along
arc-length (via `loadPaperRateDirAlongArcLen`, which rebuilds the paper's node
positions at its `dxy = 200 m` spline sampling — the source for the
interpolation). `meshGenLib.py` is not shipped here; `create.newcase` copies
it in from `scripts/meshGenLib.py`.

## Loading inputs (origin)

Per-node γ (max shear strain rate) and φ (compression angle) come from a tight
pair of files copied in from `paper.saf.A/`: `Rate_direction.txt` (rad/s, deg)
and `x1_1.txt`, `x2_1.txt`, `x3_1.txt`. These four must travel together —
`Rate_direction.txt` rows are keyed to the per-fault node positions from
spline-sampling `x*_1.txt`.

## Loading encoding (paper-faithful)

`nsmpGeoPhys.txt` per-node values:

| Col | Field | saf.gmsh.lite value |
|---|---|---|
| 6 | ftLoadMaxShear (γ) | per-node from paper Rate_direction.txt col 1 |
| 7 | ftLoadAngle (φ) | per-node from paper Rate_direction.txt col 2 |
| 8 | ftLoadWt | constant 450 |
| 9 | ftVis | per-node = `ant0·str/γ` (= paper's `ant`) |

`interstress.f90` C_mesh=3 branch computes `ant = ftVis·450/ftLoadWt =
ant0·str/γ`, matching paper.saf.A's interstress formula exactly.

## Files

```
user_defined_params.py            # case parameters
userDefinedFaultSysGeoPhys.py     # base ftPhys row template
meshgen.py                        # gmsh + paper-loading interp
user_fault_geometry_input/{ssaf,sjfn,sjfs}.gmt.txt   # SAF/SJFN/SJFS control points
README.md
```

## Resolution comparison

| | paper.saf.A | saf.gmsh.lite |
|---|---|---|
| Mesh | C_mesh=2 (Fortran structured quad) | C_mesh=3 (gmsh) |
| dxy | 200 m | 400 m |
| Fault nodes (sjfn/sjfs/ssaf) | 295 / 178 / 1769 | ~158 / ~95 / ~959 |
| Cycle 1 first interval | 471 yr | 473 yr |
| Cycle 1 magnitude | M5.74 | M5.7 |
