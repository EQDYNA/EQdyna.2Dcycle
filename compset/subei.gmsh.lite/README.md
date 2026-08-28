# subei.gmsh.lite

Coarser-resolution gmsh-meshed (C_mesh=3) Subei fault system. Mirrors
saf.gmsh.lite: spline-tangent at every gmsh fault node, paper-Fig.2 style
per-node loading plot at meshgen time. Standard build/run flow: repo root
`README.md` § Run (`create.newcase` → `python3 meshgen.py` → `case.setup` →
`bash run.sh`).

`meshgen.py` reads .gmt control points, builds a natural cubic spline y(x) per
fault (matches `meshgen1.f90` / saf.gmsh.lite; for atf, atf1+atf2 points are
merged into one x-monotonic spline), meshes with gmsh + multi-surface
decomposition for the connected junctions (`dx=0.5 km`), and for every gmsh
fault node sets `ftLoadAngle = uniformXLoadingAngle(tx, ty)` = -atan2(ty, tx)
deg. Writes
`fem_mesh_output/{vert,fac,nsmp,nsmpGeoPhys,nsmpTanLen,meshGeneralInfo}.txt`
and a loading-inputs sanity plot to `aPlots/loading_inputs.png`.
`meshGenLib.py` is not shipped here; `create.newcase` copies it in from
`scripts/meshGenLib.py`.

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
stepover), connected at junctions (atf↔dxs, atf↔sbt, sbt↔dxs); the
multi-surface gmsh decomposition handles the connected topology.

| Logical fault | gmsh trace(s) | Type | Notes |
|---|---|---|---|
| atf | atf1 + atf2 (with bridge) | strike-slip | Altyn Tagh Fault |
| dxs | dxs | strike-slip | Danghe Nan Shan |
| sbt | sbt | thrust (dip 30°) | Subei thrust |

## Files

```
user_defined_params.py, userDefinedFaultSysGeoPhys.py, meshgen.py   # case params, loading dispatch, gmsh gen
user_fault_geometry_input/{atf1,atf2,dxs,sbt}.gmt.txt   # atf segments (stepover), dxs, sbt
README.md
```

## Resolution

| | subei.gmsh.lite |
|---|---|
| Mesh | C_mesh=3 (gmsh, multi-surface) |
| dx (gmsh) | 0.5 km |
| Fault nodes (atf/dxs/sbt) | ~205 / ~99 / ~27 |
