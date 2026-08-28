# xianshuihe.gmsh.lite

Xianshuihe fault system (鲜水河断裂), eastern Tibet. **Work in progress —
geometry input only; not yet a runnable compset.**

## What is here

```
user_fault_geometry_input/
  xianshuihe_fault_trace.kml   # 1:1,000,000 Chinese active-fault database trace
plot_fault_trace.py            # parse + plot the KML, assign fault ids
map_step_overs.py              # measure and classify step-overs
xianshuihe_fault_trace.png     # output of plot_fault_trace.py
xianshuihe_step_overs.png      # output of map_step_overs.py
```

```bash
python3 plot_fault_trace.py    # figure + per-polyline table + fault id table
python3 map_step_overs.py      # figure + junction and strand step-over tables
```

Both write a 300 dpi PNG and a vector PDF, sized to a 190 mm
double-column text width.

## Fault trace

The KML holds **9 polylines across 3 Placemarks**, 124 points, 482 km of
polyline, spanning 102.44 E / 28.93 N (SE) to 100.18 E / 31.74 N (NW).
It is not one continuous line — the trace splits into sub-parallel strands
between about 30.0 and 30.5 N, the Yalahe / Selaha / Zheduotang strands
near Kangding.

**One fault per polyline, with a single deliberate merge.** The polylines
*are* the segmentation: merging them wholesale into one through-going
trace would erase the step-overs that Qiao et al. use as their maturity
metric, and would force invented geometry across the junction gaps. Ids
are assigned by `assign_fault_ids()` from `MERGE_GROUPS`, running SE to NW
— note the rotated x axis increases toward the NW, so the ascending-x sort
puts the southeastern-most group first.

**8 faults:**

| fault | polylines | length | fault | polylines | length |
|---|---|---|---|---|---|
| ft1 | seg1 + seg2 | 74 km | ft5 | seg5 | 71 km |
| ft2 | seg0 | 43 km | ft6 | seg6 | 54 km |
| ft3 | seg3 | 27 km | ft7 | seg7 | 75 km |
| ft4 | seg4 | 75 km | ft8 | seg8 | 63 km |

### The seg1 + seg2 merge

`MERGE_GROUPS` merges seg1 into seg2, and nothing else. Two reasons, both
needed — either alone would not justify it:

1. **They are not separated by a step-over.** Their endpoints meet at
   **0.26 km**, the tightest junction in the whole dataset and far below
   the 3.8 km step-over cut-off. This is a break in the digitisation, not
   a structure.
2. **seg1 alone is too short to model.** At 15 km and the `.gmsh.lite`
   resolution of dxy = 400 m it carries ~38 fault nodes, against a 2 km
   nucleation-patch radius in the SAF setup — near the floor of what the
   code resolves meaningfully. Merged, ft1 is 74 km.

No other junction qualifies: the next-tightest are 0.00 km (seg0→seg4)
and 1.62 km (seg5→seg6), but both connect polylines that are already long
enough to stand alone, and merging across them would span the splay zone.

`THROUGH_GOING_ORDER` in `plot_fault_trace.py` is **not** a fault. It is
only the order in which the polylines succeed one another along the trace,
which `map_step_overs.py` needs to know which pairs are neighbours. It is
unaffected by the merge.

## Step-overs

Qiao et al. use step-over density as a fault-maturity proxy, counting
step-overs wider than 1% of fault length: 4/350 = 0.011 per km for the
Xianshuihe against 1/550 = 0.002 for Yushu-Ganzi. `map_step_overs.py`
measures each junction (along-strike gap, fault-normal width against the
**local** strike, left/right sense) and the separation between the
parallel splay strands.

Sense is converted to releasing/restraining with the **sinistral**
convention: on this left-lateral fault a *left* step opens and a *right*
step closes — the opposite of the dextral case most references describe.

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
| seg3-seg5 | 20 km | 11.4 - 15.0 km | no (= the other two summed) |

**Negative result worth recording: this 1:1M compilation does not resolve
the step-overs the paper counts.** Every junction offset is 0.0-3.4 km,
all below the 3.8 km cut-off; only the splay strands clear it. That gives
2 step-overs, density 0.005 per km against the paper's 0.011.
Reproducing their metric needs a finer source (Chevalier et al. 2018 /
Xu et al. 2016, or 1:250k mapping).

Note the largest junction offset — 12.4 km at x = 41 km — is almost
entirely along-strike (2.3 km fault-normal). It is an unmapped stretch of
trace, not a step-over.

## Loading (not yet built)

Source: Qiao, Zhou & Zhang (2022), *EPSL* 596, 117799,
<https://doi.org/10.1016/j.epsl.2022.117799>. Their Fig. 4 gives, along
strike, everything the `nsmpGeoPhys.txt` convention needs except the
loading amplitude itself:

| `nsmpGeoPhys.txt` column | source |
|---|---|
| `ftType` | left-lateral strike-slip throughout → `1` (`FT_TYPE`) |
| `ftDip` | **90 for every fault** (`FT_DIP_DEG`). Qiao et al. Fig. 4e give 75–90 S over the Xianshuihe proper; the model takes all faults vertical. |
| `ftLoadMaxShear` | **not directly given** — see below |
| `ftLoadAngle` | `-999` (auto) |
| `ftLoadWt` | `450` |
| `ftVis` | not constrained; tune as for `gulang.gmsh.lite` |

Locking depth from Fig. 4d is 15.1 +/- 6.6 km along strike, against 22 km
for the SAF — that sets `seismoDepth`.

### Converting slip rate to `ftLoadMaxShear`

Qiao et al. give a fault **slip rate** (Fig. 4b: ~5–6 mm/yr on
Dangjiang–Yushu–Ganzi rising to ~12–13 mm/yr on the Xianshuihe), while
`ftLoadMaxShear` is a far-field shear **strain rate** in s^-1.

There is no clean kinematic conversion, because the loading rate here is
not a velocity boundary condition. `interstress.f90` uses it to set an
asymptotic *stress*, `rs = rd * cos(2*theta) * ant`: with the SAF's
`ant0 = 8.4e21 Pa s` and `1.427e-14 s^-1` that is ~120 MPa. The long-term
slip rate is then an emergent property of the model, not an input. Trying
to back out a domain width from `gamma = V / 2W` gives a nonsense answer
here (the SAF mesh is ~1300 km across, which would imply ~570 mm/yr).

Calibrate instead. Long-term slip rate is linear in the loading rate at
fixed viscosity and friction, so one run pins the constant:

1. First guess by scaling the SAF anchor. The SAF uses
   `ftLoadMaxShear = 1.427e-14 s^-1` and produces a peak modeled slip
   rate of ~38 mm/yr. For a ~13 mm/yr Xianshuihe target:
   `1.427e-14 * 13/38 ~= 4.9e-15 s^-1`.
2. Run a few hundred cycles and measure the modeled long-term slip rate
   with `scripts/plot_saf_figure3.py`.
3. Rescale linearly: `new = old * (target rate / modeled rate)`, and
   confirm on a second run.

The along-strike *shape* of the profile — the ~2x SE-ward increase — is
imposed per node through `ftLoadMaxShear`, which is a per-node column, so
the interesting part of Fig. 4b is reproducible directly. That
along-strike contrast is the scientific target of the case.

## Still to do

- Decide how the splay strands ft3/ft4/ft5 should interact at the mesh
  level (the `subei.gmsh.lite` multi-surface junction handling is the
  precedent).
- Choose the model extent: the full 900 km Yushu–Ganzi–Xianshuihe system,
  or the ~350 km Xianshuihe proper (Qiao et al. put the segment boundary
  near 100.7 E).
- Pick the origin and rotation. The trace's principal axis is
  **122.9 deg**, computed in `plot_fault_trace.py`; the SAF equivalent is
  the `x0=3.7e5, y0=3.8e6, theta=40` convention.
- Write `user_defined_params.py`, `userDefinedFaultSysGeoPhys.py`,
  `meshgen.py` following `gulang.gmsh.lite`.
