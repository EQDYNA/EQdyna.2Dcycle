# xianshuihe.gmsh.lite

Xianshuihe fault system (鲜水河断裂), eastern Tibet — runnable C_mesh=3 compset,
per-node loading from GSRM v2.1. Standard build/run flow: repo root
`README.md` § Run.

## Contents

```
user_fault_geometry_input/{xianshuihe_fault_trace.kml (1:1,000,000 Chinese active-fault DB), ft1..ft7.gmt.txt}, strain_rate_input/README.md
plot_fault_trace.py            # KML parse+plot; shared lib used by the 3 scripts below
map_step_overs.py, export_fault_geometry.py, fetch_strain_rate.sh, strain_rate_loading.py, apply_strain_loading.py   # see Step-overs / Loading
meshgen.py, user_defined_params.py, userDefinedFaultSysGeoPhys.py
```

Committed figures: `xianshuihe_fault_trace.png` (`plot_fault_trace.py`),
`xianshuihe_strain_map.png` (`strain_rate_loading.py`),
`xianshuihe_strain_figure2.png` (`strain_rate_loading.py --case <case_dir>`).
Everything else, including `xianshuihe_strain_loading.csv`, is
regenerable/gitignored.

## Fault trace

KML: **9 polylines / 3 Placemarks**, 124 points, 482 km, 102.44 E/28.93 N (SE)
to 100.18 E/31.74 N (NW); splits into sub-parallel Yalahe/Selaha/Zheduotang
strands near Kangding (30.0-30.5 N). Polylines *are* the segmentation —
merging wholesale would erase the step-overs Qiao et al. use as their maturity
metric and invent geometry across junctions. Ids from `assign_fault_ids()` via
`MERGE_GROUPS`, SE to NW.

**7 faults**, `MERGE_GROUPS = [[1,2], [0,4], [3], [5], [6], [7], [8]]`,
`FT_DIP_DEG = 90`, `FT_TYPE = 1` (all vertical, left-lateral):

| fault | polylines | length | fault | polylines | length |
|---|---|---|---|---|---|
| ft1 | seg1 + seg2 | 74.5 km | ft5 | seg6 | 54.5 km |
| ft2 | seg0 + seg4 | 118.1 km | ft6 | seg7 | 74.7 km |
| ft3 | seg3 | 26.6 km | ft7 | seg8 | 63.4 km |
| ft4 | seg5 | 70.9 km | | | **482.7 km** |

**The two merges** (nothing else is merged):

| merge | gap | justification |
|---|---|---|
| seg1+seg2 → ft1 (74.5 km) | 0.26 km | tightest junction, far below 3.8 km cut-off (digitisation break, not structure); *and* seg1 alone is 15 km → at dxy=400 m only ~38 nodes vs 2 km nucleation-patch radius, near resolution floor. Neither reason alone sufficient. |
| seg0+seg4 → ft2 (118.1 km) | 0.00 km | collinear across the join; ft2 = SE trunk |

`THROUGH_GOING_ORDER` is not a fault — polyline adjacency order for
`map_step_overs.py`, unaffected by merges.

## Step-overs

Qiao et al.'s maturity metric: step-overs > 1% of fault length — Xianshuihe
4/350 = 0.011/km vs Yushu-Ganzi 1/550 = 0.002/km. `map_step_overs.py` is the
provenance for the tables below. Sense → releasing/restraining via the
**sinistral** convention: *left* step opens, *right* step closes (opposite of
the dextral case most references use).

| junction | gap (km) | width (km) | sense | counted |
|---|---:|---:|---|---|
| 1 | 0.3 | 0.0 | right | no |
| 2 | -1.1 | 0.5 | right | no |
| 3 | 0.0 | 0.0 | right | no |
| 4 | 12.4 | 2.3 | right | no |
| 5 | -5.0 | 2.5 | left | no |
| 6 | -2.8 | 3.4 | left | no |

| strand pair | overlap | separation min-max | counted |
|---|---:|---:|---|
| seg3-seg4 | 26 km | 2.9 - 8.7 km | yes |
| seg4-seg5 | 55 km | 3.2 - 12.9 km | yes |
| seg3-seg5 | 20 km | 11.4 - 15.0 km | no (= other two summed) |

**Negative result: this 1:1M compilation does not resolve the paper's
step-overs.** All offsets 0.0-3.4 km, below the 3.8 km cut-off; only splay
strands clear it → 2 step-overs, density 0.005/km vs paper's 0.011/km. Needs a
finer source (Chevalier et al. 2018 / Xu et al. 2016, or 1:250k mapping).
Largest offset, 12.4 km at x=41 km, is almost entirely along-strike (2.3 km
fault-normal) — an unmapped trace gap, not a step-over.

## Mesh

`meshgen.py`: `saf.gmsh.lite` **embedded** scheme (one surface, fault lines
via `gmsh.model.mesh.embed`), not `gulang.gmsh.lite`'s connector-line pattern
— would invent geometry across the step-overs the case studies. Pure quads via
`Mesh.SubdivisionAlgorithm = 1`; overlapping strands preserved, no trimming.
Current mesh (dxy=400 m, `dxAtBoundary`=20 km): 56,774 nodes, **54,068 quads,
0 triangles**, 0 orphaned split nodes, 0 `checkMeshQuality.py` failures,
interior angles 32.0-154.6 deg, aspect ratio max 4.43, ~9 s. `python3
scripts/plotMeshFaults.py <case_dir>` — per-tip view matters, every orphaned
split node found in this project sat at a tip.

## Loading

Per-node from **GSRM v2.1** (Kreemer et al. 2014), sampled at mesh fault nodes
— must run after meshing:

```bash
bash fetch_strain_rate.sh                                       # GSRM v2.1 regional cut
python3 export_fault_geometry.py                                # polylines -> ftN.gmt.txt
# create.newcase / case.setup / meshgen.py — root README.md § Run
python3 <compset>/strain_rate_loading.py --case .                # sample MESH nodes (2581)
python3 <compset>/apply_strain_loading.py --case . --target-stress 100e6   # patch cols 6,7,9
```

Trap: `--case`-less `strain_rate_loading.py` samples the 123 KML points, not
the 2581 mesh nodes.

`apply_strain_loading.py` patches whichever copy of `nsmpGeoPhys.txt` it finds:
`fem_mesh_output/` on a freshly meshed case, the case root once `run.sh` has
copied the mesh up. Either is fine — `run.sh` keeps a working copy newer than
its `fem_mesh_output` source rather than overwriting it. Before that fix it
re-meshed on every launch and silently reverted the patch, and the run went
ahead on default uniform loading with no trace in the log.

| column | value |
|---|---|
| `ftLoadMaxShear` (6) | per-node max shear strain rate from GSRM |
| `ftLoadAngle` (7) | per-node angle of compression = local fault strike − max-shear direction |
| `ftLoadWt` (8) | 450 (unchanged) |
| `ftVis` (9) | `TARGET / gamma_i` — see below |

| fault | nodes | length | gamma (nanostrain/yr) | angle of compression (deg) |
|---|---:|---:|---:|---:|
| ft1 | 457 | 74.5 km | 39-97 | -21.4 .. +13.1 |
| ft2 | 613 | 118.1 km | 64-164 | -16.7 .. +41.5 |
| ft3 | 141 | 26.6 km | 60-124 | -9.2 .. +24.0 |
| ft4 | 371 | 70.9 km | 78-134 | -5.2 .. +57.4 |
| ft5 | 287 | 54.5 km | 96-165 | -13.0 .. +5.4 |
| ft6 | 391 | 74.7 km | 43-202 | -32.8 .. +18.5 |
| ft7 | 321 | 63.4 km | 74-117 | -20.4 .. +8.4 |

**Sign convention** (Liu et al. 2022, Fig. 2d): angle of compression = local
fault strike minus max-shear direction; `interstress.f90` applies `rn =
-gamma*sin(2*phi)*ant` — positive phi clamps, negative unclamps; reversed sign
pulls the fault apart. `strike`/`shear_dir` are azimuths in the local rotated
frame (principal axis 122.9 deg) — **rotation cancels in the difference**,
immune to frame choice. GSRM tensor **components** are interpolated, never the
principal angle (pi-periodic; interpolating it directly spikes +/-90 deg at
the wrap).

**`ftVis` and TARGET.** `ftVis_i = TARGET/gamma_i` (SAF construction, `ant =
ant0*str/rd`): constant rate x viscosity → uniform asymptotic shear stress,
gamma sets only the approach timescale. TARGET must not exceed ambient normal
stress: with `ambientnorm = -100 MPa`, `ns = -N - T*sin(2*phi)`, `strength =
fric_fs*|ns|`, `shs = T*cos(2*phi)`, tensile once `T*|sin 2phi| > N`
(negative-phi side):

| T/N | tensile nodes | strength < 5 MPa | can nucleate | max ns |
|---:|---:|---:|---:|---:|
| 0.90 | 0% | 0% | 90% | -18.0 MPa |
| 1.00 | 0% | 0.4% | 93% | -8.9 MPa |
| 1.20 | 1.5% | 2.4% | 95% | +9.3 MPa |

SAF uses T/N=1.20 because its angle of compression is coherently positive (91%
of nodes, mean +6.8 deg); Xianshuihe straddles zero (43% positive, mean -0.2
deg). **Use T = 90-100 MPa here.** Backstop: `interstress.f90` now caps
effective normal stress at `minnorm = -10 MPa` (matches `faulting.f90`'s
dynamic cap; previously `abs(ns)`, giving a tensile node strength proportional
to how tensile it was).

| compset | ftLoadMaxShear | nanostrain/yr | ftVis (Pa s) |
|---|---|---|---|
| `saf.gmsh.lite` | varies per node | 143-536 | 7.1e21-2.6e22 |
| `subei.gmsh.lite` | 1.427e-14 uniform | 450 | 8.4e21 |
| `gulang.gmsh.lite` | 1.427e-14 uniform | 450 | 5.0e21 |
| `xianshuihe` (GSRM) | varies per node | 39-202 | 1.6e22-8.0e22 |

Only SAF and this case derive loading from a strain field; `subei`/`gulang`
carry the SAF default 450 nanostrain/yr unchanged, above the SAF's own max
everywhere along their faults. Sanity check (`V=2*W*gamma`, 50 km half-width):
5.4 mm/yr (ft1), 11-12 mm/yr (Xianshuihe proper) vs Qiao et al.'s 5-6/12-13
mm/yr — independent InSAR product, not fitted (same conversion: SAF 45 vs
observed 34, order-of-magnitude only, W free). `ftVis` reaches 8.0e22 Pa s at
the lowest strain rate, above SAF's 2.6e22 max (constant stress-product
construction, not a rheological claim). Gamma doesn't gate rupture (cancels
from asymptotic stress); sets only the clock (`tau=ftVis/mu`: 15-25 kyr here
vs ~8 kyr SAF).

## Reference

Qiao, Zhou & Zhang (2022), *EPSL* 596, 117799,
<https://doi.org/10.1016/j.epsl.2022.117799>. Fig. 4d: locking depth 15.1 +/-
6.6 km along strike (vs 22 km SAF); Fig. 4e: dips 75-90 S over the Xianshuihe
proper (model rounds to vertical).

## Still to do

- Long runs for a converged catalogue; recurrence interval currently a few years (small-event dominated).
- Calibrate against Qiao et al.'s along-strike slip-rate profile with `scripts/plot_saf_figure3.py`.
- Compare the ~350 km Xianshuihe proper against the full system (segment boundary near 100.7 E, Qiao et al.).
