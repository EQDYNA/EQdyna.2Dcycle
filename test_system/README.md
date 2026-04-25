# test_system

Lightweight regression check for EQdyna.2Dcycle.

## What it does

For each compset listed in `testNameList.py`:
1. Build the binary from `src/`.
2. `create.newcase` → `meshgen.py` → `case.setup` → `run.sh`.
3. Diff `aRawSimuData/{cyclelog,interval,totalop}.txt*` against
   `reference.results/<case>/`.

Any mismatch beyond `1e-3` → non-zero exit.

## Running

```bash
# Full pipeline (compile + run + verify)
python3 -m test_system.test_all

# Just verify already-run results
python3 test_system/verify.test.py

# Quick smoke (compile + 1 short case + check totalop.txt1 exists)
python3 test_system/smoke.py
```

## Adding a new test case

1. Add the compset name to `testNameList.py`.
2. Run the case once manually, then copy
   `aRawSimuData/{cyclelog,interval,totalop}.txt*` into
   `reference.results/<case>/`.
3. Re-run `test_all` to confirm it now passes.

## Notes

- `paper.saf.A` is intentionally NOT in the regression set — it takes
  hours and is a research artifact. Use it as a manual reference run
  when validating physics changes.
- Tolerance is `1e-3` element-wise. Tighten in `verify.test.py` if you
  need stricter comparison.
