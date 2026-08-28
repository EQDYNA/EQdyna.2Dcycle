# gulang.gmsh.lite

Coarser-resolution gmsh-meshed (C_mesh=3) Gulang fault system (ft1..ft5, all
left-lateral strike-slip, vertical dip). Mirrors
saf.gmsh.lite/subei.gmsh.lite: spline-tangent at every gmsh fault node,
paper-Fig.2-style per-node loading plot at meshgen time. Standard build/run
flow: repo root `README.md` § Run (`create.newcase` → `python3 meshgen.py` →
`case.setup` → `bash run.sh`). `meshgen.py` reads .gmt control points, builds
a natural cubic spline y(x) per fault (matches `meshgen1.f90`), meshes with
gmsh (`dx=0.3 km`), and writes per-node loading via `uniformXLoadingAngle` in
`scripts/meshGenLib.py` (copied in by `create.newcase`, not shipped here),
writes
`fem_mesh_output/{vert,fac,nsmp,nsmpGeoPhys,nsmpTanLen,meshGeneralInfo}.txt`,
and saves a loading-inputs sanity plot to `aPlots/loading_inputs.png`.

## Loading model (default: uniform max-shear strain rate along +x)

Gulang has no paper-A-style `Rate_direction.txt` source, so per-node loading
is derived from fault geometry alone:

```
ftLoadAngle (deg) = atan2(ty, tx)         # = fault_strike - 0 (max-shear along +x)
ftLoadMaxShear   = 1.427e-14 s^-1         # uniform reference rate γ
ftVis            = 5.0e21 Pa·s            # see "ftVis tuning" below
ftLoadWt         = 450
ftType           = 1 (left-lateral strike-slip), ftDip = 90 (vertical)
```

The 45° upper-clamp in `interstress.f90` keeps vertical/perpendicular faults
passive (`cos(2·45°)=0` → no shear loading, no nucleation).

### ftVis tuning (5e21 instead of paper's 8.4e21)

Paper Model A uses ftVis = 8.4e21 Pa·s, giving asymptotic shear amplitude
`γ·ant ≈ 120 MPa`. For paper.A this is fine because φ is mostly 0–30° →
`sin(2φ)` small → normal-stress changes small → fault stays compressive. For
gulang, φ spans ±45° (after clamping) → `|sin(2φ)|` reaches 1 → normal-stress
changes can swing ±120 MPa around the −100 MPa ambient → tensile normal stress
→ friction breaks → NaN.

| ftVis | outcome |
|---|---|
| 8.4e21 | tensile normal at large strikes → NaN at fault 3 |
| 1.0e21 | only 14 MPa shear amplitude, never reaches the 50 MPa nucleation threshold |
| **5.0e21** (current) | rs∞ ≈ 71 MPa (above threshold), \|rn∞\| ≈ 71 MPa (below ambient) — runs, 10 cycles in test |

### ft3 / ft4 trim (geometric overlap fix)

The original `ft3.gmt.txt` ended at `(10.09, -0.45) km` while `ft4.gmt.txt`
starts at `(9.21, -0.38) km` — they overlap ~0.9 km in x with only ~70 m
vertical separation, well below `dx = 0.3 km`. gmsh produces sliver elements
between them, giving some fault nodes zero mass → NaN at the first
dynamic-rupture step. Fix: `ft3.gmt.txt` was trimmed to end at `(8.96, -0.36)
km` (drop the last 3 control points), leaving a ~0.25 km clean gap before ft4;
the original tail is preserved in git history.

## Files

```
user_defined_params.py              # case parameters
userDefinedFaultSysGeoPhys.py       # gulang loading dispatch (uniform per fault)
meshgen.py                          # gmsh + spline-tangent + uniform-x loading
user_fault_geometry_input/ft1.gmt.txt .. ft5.gmt.txt   # fault control points
README.md
```

## Resolution

| | gulang.gmsh.lite |
|---|---|
| Mesh | C_mesh=3 (gmsh) |
| dx (gmsh) | 0.3 km |
| Fault nodes (ft1..ft5) | 205 / 66 / 111 / 185 / 125 |
