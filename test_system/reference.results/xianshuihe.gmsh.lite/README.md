# Frozen reference — xianshuihe.gmsh.lite, 5 cycles

A solver regression fixture. It pins the **dynamic rupture and interseismic
solver**, not the mesher or the loading pipeline: the mesh and the per-node
loading are shipped here as frozen inputs, so a gmsh version change or a GSRM
re-download cannot move the answer. Meshing and loading are covered separately
by `scripts/checkMeshQuality.py` and the compset's own scripts.

## Provenance

| | |
|---|---|
| Compset | `xianshuihe.gmsh.lite`, 7 faults, 2581 fault nodes |
| Mesh | dxy = 400 m, pure quad, 54,068 elements, 0 orphans |
| Loading | GSRM v2.1, per-node, `--target-stress 100e6` (T/N = 1.0) |
| Binary | `run_eqdyna2d_2.1.1` |
| Threads | **`OMP_NUM_THREADS=1`** — see R30; the answer is thread-count dependent |
| Cycles | 1-5 |
| Intervals | 1486, 53, 38, 180, 117 yr |

## Files

| File | Content |
|---|---|
| `inputs.tar.gz` | complete frozen input state: mesh, `nsmpGeoPhys.txt`, `FE_*.txt`, `user_defined_params.py` |
| `totalop.txt1.gz` | expected output, 5 cycles x 2581 nodes x 5 columns |
| `interval.txt1`, `cyclelog.txt1` | expected recurrence intervals and cycle log |
| `md5sums.txt` | checksums for the two archives |

## Why this case

It is the only compset exercising more than four faults, which is where three
silent Fortran defects lived until v2.1.0 — the `faulting.f90` nucleation
lookup with no else clause past fault 4, the `interstress.f90` fault-tip and
nucleation loops hardcoded to three faults, and uncapped tensile strength.
None of them could be caught by a 3-fault SAF case. Nucleation in these five
cycles lands on ft5, which is exactly the range the old code could not
resolve.

## Regenerating

Only when a change is *intended* to alter results. Record why in the commit.

```bash
cd work/xsh_ref_gen && OMP_NUM_THREADS=1 ./run_eqdyna2d_<version>
# then re-pack inputs.tar.gz, totalop.txt1.gz, interval.txt1, cyclelog.txt1, md5sums.txt
```
