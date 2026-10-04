# Testing and reproducibility

## The test suite

```bash
python3 test_system/test_conventions.py     # fast checks on code and docs conventions
python3 test_system/smoke.py                # compile and run one cycle
python3 test_system/verify_xianshuihe.py    # solver regression: 5 cycles, bit-exact
python3 test_system/verify_stf.py           # rupture-time-history output
python3 -m test_system.test_all             # full pipeline
```

`verify_xianshuihe.py` runs five cycles of a 7-fault model from frozen inputs
and requires the output to match the stored reference bit for bit. The mesh
and loading are frozen too, so the test checks the solver alone.

`verify_stf.py` checks that turning the time-history output on leaves the
solution unchanged, and that the time histories agree with the end-of-event
output, integrate to the final slip, include every node that ruptured, carry
the right event times, and reproduce the catalogue moment.

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
