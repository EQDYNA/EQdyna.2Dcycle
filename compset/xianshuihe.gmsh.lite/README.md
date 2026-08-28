# xianshuihe.gmsh.lite

Xianshuihe fault system (鲜水河断裂), eastern Tibet. **Work in progress —
geometry input only; not yet a runnable compset.**

## What is here

```
user_fault_geometry_input/
  xianshuihe_fault_trace.kml   # 1:1,000,000 Chinese active-fault database trace
plot_fault_trace.py            # parse + plot the KML, propose the fault decomposition
xianshuihe_fault_trace.png     # output of the above
```

Run `python3 plot_fault_trace.py` to regenerate the figure and print the
per-polyline table and the gap report.

## Fault trace

The KML holds **9 polylines across 3 Placemarks**, 124 points, 482 km of
polyline, spanning 102.44 E / 28.93 N (SE) to 100.18 E / 31.74 N (NW).
It is not one continuous line — the trace splits into sub-parallel strands
between about 30.0 and 30.5 N, the Yalahe / Selaha / Zheduotang strands
near Kangding.

Proposed decomposition (`FAULT_DECOMPOSITION` in the script):

| fault | polylines | length | gaps |
|---|---|---|---|
| ft1 through-going | 1, 2, 0, 4, 6, 7, 8 | 385 km | 12.6 km at x=41, 5.6 at x=108, 4.4 at x=177 |
| ft2 parallel strand | 5 | 71 km | none |
| ft3 short splay | 3 | 27 km | none |

**The gaps are unresolved and have to be closed when the mesh is built.**
Neither strand gives a gap-free path through the Kangding splay zone:
seg4 joins seg0 exactly (0.00 km) but stops 12.6 km short of seg6, while
seg5 reaches seg6 (1.6 km) but starts 19.8 km off seg0. ft1 takes seg4 as
the smaller-maximum-gap option. The script draws each gap dotted with an
x rather than bridging it, so the figure never implies a continuity the
data does not have.

Topologically this is the same shape as the SAF compset — one
through-going fault plus splay strands — so `subei.gmsh.lite`'s
multi-surface junction handling is the relevant precedent.

## Loading (not yet built)

Source: Qiao, Zhou & Zhang (2022), *EPSL* 596, 117799,
<https://doi.org/10.1016/j.epsl.2022.117799>. Their Fig. 4 gives, along
strike, everything the `nsmpGeoPhys.txt` convention needs except the
loading amplitude itself:

| `nsmpGeoPhys.txt` column | source |
|---|---|
| `ftType` | left-lateral strike-slip throughout → `1` |
| `ftDip` | Fig. 4e: 60–90 S (Dangjiang/Yushu), 70–90 N (Ganzi), 75–90 S (Xianshuihe) |
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

- Close the three ft1 gaps and settle the splay-zone junction.
- Choose the model extent: the full 900 km Yushu–Ganzi–Xianshuihe system,
  or the ~350 km Xianshuihe proper (Qiao et al. put the segment boundary
  near 100.7 E).
- Pick the origin and rotation. The trace's principal axis is
  **122.9 deg**, computed in `plot_fault_trace.py`; the SAF equivalent is
  the `x0=3.7e5, y0=3.8e6, theta=40` convention.
- Write `user_defined_params.py`, `userDefinedFaultSysGeoPhys.py`,
  `meshgen.py` following `gulang.gmsh.lite`.
