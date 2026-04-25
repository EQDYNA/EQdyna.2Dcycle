# gulang.gmsh.lite

Coarser-resolution gmsh-meshed (C_mesh=3) Gulang fault system. Mirrors the
saf.gmsh.lite + subei.gmsh.lite framework: spline-tangent at every gmsh
fault node, paper-Fig.2-style per-node loading plot at meshgen time.

## Loading model (default: uniform max-shear strain rate along +x)

Gulang has no paper-A-style `Rate_direction.txt` source, so per-node
loading is derived from fault geometry alone using `uniformXLoadingAngle`
in `scripts/meshGenLib.py`:

```
ftLoadAngle (deg) = atan2(ty, tx)         # = fault_strike - 0 (max-shear along +x)
ftLoadMaxShear   = 1.427e-14 s^-1         # uniform reference rate γ
ftVis            = 5.0e21 Pa·s            # see "ftVis tuning" below
ftLoadWt         = 450
ftType           = 1 (left-lateral strike-slip), ftDip = 90 (vertical)
```

The 45° upper-clamp in `interstress.f90` keeps vertical/perpendicular
faults passive (`cos(2·45°)=0` → no shear loading, no nucleation).

### ftVis tuning (5e21 instead of paper's 8.4e21)

Paper Model A uses ftVis = 8.4e21 Pa·s, giving asymptotic shear
amplitude `γ·ant ≈ 120 MPa`. For paper.A this is fine because the φ
values are mostly 0–30° → `sin(2φ)` is small → normal-stress changes
are small → fault stays compressive.

For gulang the per-node φ spans ±45° (after clamping) → `|sin(2φ)|`
reaches 1 → normal-stress changes can swing ±120 MPa around the
−100 MPa ambient → tensile normal stress → friction breaks → NaN.

Empirically:
- `ftVis = 8.4e21`: tensile normal at large strikes → NaN at fault 3.
- `ftVis = 1.0e21`: only 14 MPa shear amplitude, never reaches the
  50 MPa nucleation threshold.
- **`ftVis = 5.0e21`** (current): rs∞ ≈ 71 MPa (above threshold) and
  |rn∞| ≈ 71 MPa (below ambient) — simulation runs, 10 cycles in test.

### ft3 / ft4 trim (geometric overlap fix)

The original `ft3.gmt.txt` ended at `(10.09, -0.45) km` while
`ft4.gmt.txt` starts at `(9.21, -0.38) km` — they overlap ~0.9 km in
x with only ~70 m vertical separation, well below `dx = 0.3 km`.
gmsh produces sliver elements between them → some fault nodes get
zero mass → NaN at the first dynamic-rupture step.

Fix: `ft3.gmt.txt` was trimmed to end at `(8.96, -0.36) km` (drop the
last 3 control points), leaving a ~0.25 km clean gap before ft4. The
original tail is preserved in git history if needed.

## Fault system

Five embedded fault traces, all left-lateral strike-slip with vertical dip
(ft1..ft5 from `user_fault_geometry_input/`).

## Mesh + loading flow

Same as saf.gmsh.lite / subei.gmsh.lite:
1. `meshgen.py` reads .gmt control points.
2. Builds natural cubic spline y(x) per fault (matches `meshgen1.f90`).
3. gmsh meshes 2D domain with embedded fault polylines; element size `dx=0.3 km`.
4. For every gmsh fault node:
   - tangent = analytical spline derivative (smooth, no chord kinks).
   - ftLoadAngle = `uniformXLoadingAngle(tx, ty)` = -atan2(ty, tx) deg.
   - other ftPhys fields from `userDefinedFaultSysGeoPhys.defineSysPhys`.
5. Writes `fem_mesh_output/{vert,fac,nsmp,nsmpGeoPhys,nsmpTanLen,meshGeneralInfo}.txt`.
6. Saves loading-inputs sanity-check plot to `aPlots/loading_inputs.png`.

## Files

```
user_defined_params.py              # case parameters
userDefinedFaultSysGeoPhys.py       # gulang loading dispatch (uniform per fault)
meshgen.py                          # gmsh + spline-tangent + uniform-x loading
user_fault_geometry_input/
  ft1.gmt.txt .. ft5.gmt.txt        # fault control points
README.md
```

`meshGenLib.py` is **not** in this dir — copied in by `create.newcase`
from the canonical `scripts/meshGenLib.py`.

## Usage

```bash
create.newcase --work_dir work/gulang --compset gulang.gmsh.lite
cd work/gulang
python3 meshgen.py
python3 case.setup
bash run.sh
```

## Resolution

| | gulang.gmsh.lite |
|---|---|
| Mesh | C_mesh=3 (gmsh) |
| dx (gmsh) | 0.3 km |
| Fault nodes (ft1..ft5) | 205 / 66 / 111 / 185 / 125 |
