# subei.gmsh.lite

Coarser-resolution gmsh-meshed (C_mesh=3) Subei fault system. Mirrors the
saf.gmsh.lite framework: spline-tangent at every gmsh fault node, paper-Fig.2
style per-node loading plot at meshgen time.

## Loading model (default: uniform compression along +x)

Subei has no paper-A-style `Rate_direction.txt` source, so the per-node
loading angle is derived from fault geometry alone:

```
ftLoadAngle (deg) = -atan2(ty, tx)        # = (loading dir 0°) - (fault tangent)
ftLoadMaxShear   = 1.427e-14 s^-1         # uniform reference rate
ftVis            = 8.4e21 Pa·s            # paper's ant0
ftLoadWt         = 450  (default; per-fault overrides below)
```

Per-fault overrides (in `userDefinedFaultSysGeoPhys.py`):

| Fault | ftType | ftDip | ftLoadWt | Note |
|---|---|---|---|---|
| atf (x<0) | 1 (right-strike) | 90 | 450 | full loading |
| atf (x>0) | 1 | 90 | 360 (=0.8·450) | 80% loading |
| dxs | 1 | 90 | 45 (=0.1·450) | 10% loading |
| sbt | 2 (thrust) | 30 | 450 | full loading |

## Fault system

Four embedded fault traces (atf is split into atf1+atf2 across a ~2.2 km
stepover), connected at junctions (atf↔dxs, atf↔sbt, sbt↔dxs). The
multi-surface gmsh decomposition handles the connected topology.

| Logical fault | gmsh trace(s) | Type | Notes |
|---|---|---|---|
| atf | atf1 + atf2 (with bridge) | strike-slip | Altyn Tagh Fault |
| dxs | dxs | strike-slip | Danghe Nan Shan |
| sbt | sbt | thrust (dip 30°) | Subei thrust |

## Mesh + loading flow

1. `meshgen.py` reads .gmt control points.
2. Builds natural cubic spline y(x) per fault (matches `meshgen1.f90` /
   saf.gmsh.lite). For atf, atf1+atf2 control points are merged into one
   x-monotonic spline.
3. gmsh meshes 2D domain with embedded fault polylines + multi-surface
   decomposition for connected junctions; element size `dx=0.5 km`.
4. For every gmsh fault node:
   - tangent = analytical spline derivative (smooth, no chord kinks).
   - ftLoadAngle = `uniformXLoadingAngle(tx, ty)` = -atan2(ty, tx) deg.
   - other ftPhys fields from `userDefinedFaultSysGeoPhys.defineSysPhys`.
5. Writes `fem_mesh_output/{vert,fac,nsmp,nsmpGeoPhys,nsmpTanLen,meshGeneralInfo}.txt`.
6. Saves loading-inputs sanity-check plot to `aPlots/loading_inputs.png`.

## Files

```
user_defined_params.py              # case parameters
userDefinedFaultSysGeoPhys.py       # subei loading dispatch (uniform per fault)
meshgen.py                          # gmsh + spline-tangent + uniform-x loading
user_fault_geometry_input/
  atf1.gmt.txt, atf2.gmt.txt        # atf segments (stepover)
  dxs.gmt.txt
  sbt.gmt.txt
README.md
```

`meshGenLib.py` is **not** in this dir — copied in by `create.newcase`
from the canonical `scripts/meshGenLib.py`.

## Usage

```bash
create.newcase --work_dir work/subei --compset subei.gmsh.lite
cd work/subei
python3 meshgen.py
python3 case.setup
bash run.sh
```

## Resolution

| | subei.gmsh.lite |
|---|---|
| Mesh | C_mesh=3 (gmsh, multi-surface) |
| dx (gmsh) | 0.5 km |
| Fault nodes (atf/dxs/sbt) | ~205 / ~99 / ~27 |
