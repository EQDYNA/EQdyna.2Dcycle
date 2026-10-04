# Source time functions

The standard outputs keep only each event's end state. For ground-motion
simulation you need the time history, so the solver can also write, for
**every** event, the full dynamic time series at every fault node that
slipped.

## Turning it on

```python
par.outputSTF = 1      # in user_defined_params.py
par.stfEvery = 1       # keep every Nth solver step (default 1 = every 0.01 s)
par.stfVmin = 1.e-3    # m/s; write only nodes whose peak slip rate exceeds this
```

then `python3 case.setup` and run as usual. The output is one file per run
segment, `stf.bin<icstart>`. Turning it on does not change the solution:
`totalop.txt` is bit-identical with it on or off.

## What is in it

| per file | per event | per node | per node, per time step |
|---|---|---|---|
| fault-node table: fault, x, y, unit tangent, tributary length | event id (matches `catalog.csv`), time in the sequence, preceding interval, nucleation node, sample interval, duration | rupture time, peak slip rate | slip rate, signed slip rate, slip, signed slip, shear stress, normal stress, friction |

Signed quantities are along the node's fault tangent. Values are float32.

Events are not a fixed length: an event stops once the peak slip rate falls
below 1 mm/s (after at least 5 s, at most `term`), and each record carries its
own duration. Size scales with event size: a 4000-cycle San Andreas run at full
resolution is about 45 GB.

## Reading it

```bash
python3 stf_read.py . --list                       # every event, with M0 and Mw
python3 stf_read.py . --event 24 --netcdf ev24.nc  # one event to NetCDF (compressed)
python3 stf_read.py . --event 24 --npz ev24.npz    # or numpy
python3 plot_stf.py . 24                           # figure
```

```python
from stf_read import STFCase
case = STFCase("work/my_case")
ev = case.event(24)                    # dict of numpy arrays
ev["time"], ev["x"], ev["slip_rate"]   # (nt,), (nodes,), (nodes, nt)
t, mdot = case.moment_rate(24)         # moment-rate function, N m/s
```

The binary layout is documented at the top of `src/stf_output.f90` for
readers in other languages. The plot puts every event of a case on the same
axes — the whole fault system and a fixed amplitude scale — so events can be
compared directly.

## Geographic coordinates (San Andreas compsets)

The model works in a local rotated frame. For `paper.saf.A`, the frame is SCEC
CFM 5.2 geometry in UTM zone 11, shifted and rotated by 40°.
`saf_model_to_geo.py` inverts that exactly:

```bash
python3 saf_model_to_geo.py . -o fault_nodes_geo.csv   # UTM 11, lon/lat, strike, length
python3 saf_model_to_geo.py --check                    # verify against the published numbers
```

Use the UTM columns against SCEC CFM and CVM. The lon/lat assume WGS84.

## Using it as a 3D kinematic source

Each node supplies a subfault's location, strike, slip-rate function and
rupture onset; faults are vertical and strike-slip. What the 2D model does
not define is variation with depth. A common choice is to repeat each node's
slip-rate function over the seismogenic width (22 km, as the moments assume),
tapered near the surface and the base. That is a modelling decision for the
ground-motion study.

Behind the rupture front the slip rate rings with repeated short pulses; part
of that high-frequency content is likely numerical. Low-pass filter before
using it for ground motion.
