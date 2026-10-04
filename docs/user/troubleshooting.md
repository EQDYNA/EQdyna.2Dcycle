# Troubleshooting

| symptom | likely cause | what to do |
|---|---|---|
| A parameter change has no effect | the line lacks the `par.` prefix, or `case.setup` was not re-run | write `par.<name> = …`, then `python3 case.setup` |
| Every interval after the first is 1 year, with no real slip | a node is permanently above failure — often tensile normal stress from too large a loading stress | lower the target stress so T ≤ \|σ_ambient\| where the angle of compression is negative ([Model inputs](inputs.md)) |
| No earthquakes for tens of thousands of years | asymptotic shear stress below strength | raise the loading stress or check the angle of compression |
| `ERROR: NaN in slip/sliprate` | the solution diverged at the node it reports | check the mesh at that location with `checkMeshQuality.py` and `plotMeshFaults.py` |
| `checkMeshQuality.py` reports orphaned split nodes | an element on one side of a fault is missing, usually at a fault tip | re-mesh; inspect the tips with `plotMeshFaults.py` |
| Loading is uniform although you patched `nsmpGeoPhys.txt` | an older `run.sh` re-meshed and overwrote the patch | regenerate `run.sh` with `case.setup`; current versions keep the newer file |
| Results differ from a stored reference after a few cycles | different `OMP_NUM_THREADS` | rerun at the reference's thread count ([Testing](testing.md)) |
| A plot reports no events | window past the last event, or a slip threshold above every event | pick a window inside the event range, or lower `--threshold` |
| The catalogue covers only part of a restarted run | the live directory was read directly | plot from `snapshot_case.py`'s copy, which merges the segments |
| Event times restart from zero after a restart | `seqtime.txt` was missing | keep `seqtime.txt` in the case directory when restarting |
