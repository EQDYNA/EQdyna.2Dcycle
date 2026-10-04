# Running a case

## The workflow

```bash
create.newcase --work_dir work/my_case --compset <compset>
cd work/my_case
# edit user_defined_params.py
python3 case.setup        # writes FE_*.txt and run.sh from the parameters
python3 meshgen.py        # C_mesh = 3 compsets only: builds the gmsh mesh
bash run.sh               # runs in the background
```

`create.newcase --list` lists the compsets. `create.newcase` refuses to
overwrite an existing directory; add `--force` to wipe and recreate it.

Parameters flow `defaultParameters.py` → `user_defined_params.py` →
`case.setup` → `FE_*.txt`. Set a parameter in `user_defined_params.py` as
`par.<name> = <value>`; a line without the `par.` prefix is silently ignored.

## Meshing

* **C_mesh = 2** — the solver builds a structured quadrilateral mesh itself
  from the fault geometry files `x*_1.txt`. No separate step.
* **C_mesh = 3** — `meshgen.py` builds an unstructured quadrilateral mesh with
  gmsh into `fem_mesh_output/`. Check it before a long run:

```bash
python3 checkMeshQuality.py .     # exits 1 on a real defect
python3 plotMeshFaults.py .       # every fault end to end, and every fault tip
```

`checkMeshQuality.py` fails on orphaned split nodes, triangles in a quad mesh,
interior angles below 20° or above 160°, aspect ratio above 10, and
degenerate cells. Inspect the fault tips in particular — that is where an
embedded fault's mesh goes wrong.

`run.sh` does not re-mesh when a mesh already exists, and never overwrites a
working copy newer than its `fem_mesh_output` source. Set `FORCE_MESH=1` to
re-mesh on purpose.

## Threads

`run.sh` uses `OMP_NUM_THREADS` from your environment, default 1:

```bash
OMP_NUM_THREADS=4 bash run.sh
```

Results depend on the thread count — see
[Testing and reproducibility](testing.md). Use the same count as any run you
intend to compare against.

## Restarting

Every cycle refreshes `binaryop`, the full restart state, and `seqtime.txt`,
the time in the sequence. To extend a run that ended at cycle N:

```bash
sed -i 's/par\.icstart, par\.icend = .*/par.icstart, par.icend = N+1, M/' user_defined_params.py
python3 case.setup
bash run.sh
```

For C_mesh = 2 cases, `run.sh` moves `binaryop` into `aRawSimuData/` when a
run finishes; copy it back into the case directory first. The restarted
segment writes `totalop.txt<N+1>` and friends alongside the originals; the
post-processing scripts stitch the segments together.

## Monitoring a running case

```bash
tail -f run_*.log                  # solver progress, every 3 simulated seconds
wc -l interval.txt1                # cycles completed
python3 snapshot_case.py .         # consistent copy for plotting while it runs
monitor_runs.sh -i 1200 work/my_case   # re-plot every 20 minutes
```

Plot from a snapshot, not the live directory: the live files end mid-cycle,
and a restarted case keeps its segments in separate files.
