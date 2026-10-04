# Output files

When a run finishes, `run.sh` moves the raw results into `aRawSimuData/`.
`<N>` is the first cycle of the run segment (`icstart`), so a fresh run
writes `totalop.txt1` and a restart at cycle 1257 writes `totalop.txt1257`.

| file | content |
|---|---|
| `totalop.txt<N>` | per cycle, per fault node, the end-of-event state (below) |
| `interval.txt<N>` | one line per cycle: the interseismic interval before the event, yr |
| `cyclelog.txt<N>` | first and last cycle of the segment |
| `binaryop` | full restart state after the last cycle |
| `seqtime.txt` | time in the earthquake sequence after the last cycle, yr |
| `stf.bin<N>` | optional rupture time histories — see [Source time functions](source-time-functions.md) |
| `run_*.log` | solver log |

The run's inputs (`vert.txt`, `fac.txt`, `nsmp.txt`, `nsmpGeoPhys.txt`,
`meshGeneralInfo.txt`, `FE_*.txt`, …) stay alongside.

## `totalop.txt`

All cycles stacked: `ncycles × totftnode` rows, faults in order, five
columns per node:

| col | quantity | unit |
|---|---|---|
| 1 | shear stress after the event | Pa |
| 2 | normal stress after the event (negative = compression) | Pa |
| 3 | slip in the event | m |
| 4 | slip rate when the event stopped | m/s |
| 5 | rupture time (1000 = the node did not rupture) | s |

`meshGeneralInfo.txt` gives the node count per fault, so cycle k, fault j
starts at row `(k-1)·totftnode + sum(nfnode[:j])`.

## Event times

The time of event k in the sequence is the sum of the first k intervals in
`interval.txt` (the first interval runs from the start of the model).

## Catalogue

`CATALOG=1 plotRuptureDynamics` writes `aPlots/catalog.csv`, one row per
event:

| column | meaning |
|---|---|
| `eqId` | cycle number |
| `magnitude`, `moment_Nm` | Mw and seismic moment (μ = ρ·vs², 22 km seismogenic width) |
| `nuc_x_km`, `nuc_y_km`, `nuc_ft` | nucleation point and fault (0-based) |
| `rup_dur_s` | time for the rupture front to finish |
| `peak_slip_m` | largest slip in the event |

Magnitudes use the model's own rigidity, 3.20 × 10¹⁰ Pa. Liu et al. (2022)
Figures 6 and 9 use 3.675 × 10¹⁰ Pa, so their magnitudes sit 0.04 Mw higher.
