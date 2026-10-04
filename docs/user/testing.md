# Testing and reproducibility

## The test suite

One runner, four tiers:

```bash
python3 test_system/run.py               # unit + regression (seconds)
python3 test_system/run.py unit          # pytest tests of the Python tools
python3 test_system/run.py regression    # convention checks, one per past failure
python3 test_system/run.py smoke         # build, mesh and run one cycle (minutes)
python3 test_system/run.py e2e           # bit-exact solver checks against frozen references
python3 test_system/run.py all           # everything; run before a release
```

| tier | what it checks | time |
|---|---|---|
| unit | the readers, the geographic conversion, snapshot merging, the parameter reference, the published-archive reader | < 1 s |
| regression | code and documentation conventions that once failed silently | ~3 s |
| smoke | the Fortran builds, a gmsh case meshes, one cycle runs and writes output | a few minutes |
| e2e | five cycles of a 7-fault model match the stored reference bit for bit; the rupture-time-history output leaves the solution unchanged and agrees with it | a few minutes |

Continuous integration runs unit, regression and smoke on every push and pull
request. The e2e references were made on the project's own machine; the build
uses CPU-specific instructions, so on other hardware e2e can differ in the
last bit without anything being wrong. Run it locally before a release.

## Reproducibility depends on the thread count

Runs are reproducible only at a **fixed `OMP_NUM_THREADS`**. Parts of the
solver accumulate forces in parallel, so the order of additions — and the last
bit of the result — depends on thread scheduling. An earthquake sequence is
chaotic: that last-bit difference grows into a visibly different catalogue
within a few cycles.

At one thread, both San Andreas compsets reproduce their stored references
exactly. A run at a different thread count matches for three or four cycles
and then diverges; that is expected, not a regression. Record the thread count
with any result you may want to reproduce.
