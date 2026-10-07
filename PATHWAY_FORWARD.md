# Pathway forward — EQdyna.2Dcycle

Present-tense status board (invariant 12 / R37): every open issue, to-do and
standing claim, one row, with a priority the work is taken in. State is not
priority — re-prioritise by editing the Priority column only. Blank
"Last checked" = never audited, stays blank. A claim with no command is
remembered, not verified, and is marked so.

| Item | Priority | State | Recheck | Last checked | Command | Route |
|---|---|---|---|---|---|---|
| Release v2.3.0: CHANGELOG `[Unreleased]` holds the `slipSense`/ftType fix; `VERSION` says 2.2.0 but that version was never tagged (latest tag `v2.1.1`) — fold the untagged 2.2.0 content into 2.3.0 rather than tagging 2.2.0 after the fact | P1 | open | per release | 2026-10-07 | `git describe --tags --abbrev=0` (→ `v2.1.1`) vs `cat VERSION` (→ `2.2.0`) | haruto-nakamura |
| `gulang.gmsh.lite` was left-lateral but simulated right-lateral (solver never read `ftType`) — fixed with `par.slipSense = -1` | P1 | MEASURED (PR #5, `3caff46`): 233 cycles, shear stress 100% negative (96,744/96,744 rows, mean -48.5 MPa), reproduced independently from `totalop.txt1` | per release | 2026-10-07 | `grep -n slipSense compset/gulang.gmsh.lite/user_defined_params.py` (→ `par.slipSense = -1`) | lars-eriksson |
| `xianshuihe.gmsh.lite` was left-lateral but simulated right-lateral; fixed with `par.slipSense = -1` and a rederived TARGET=90 MPa (the prior 90-100 MPa note was built for the opposite sign convention and stalled at simulated year ~2.27M with zero ruptures before the redo) | P1 | MEASURED (PR #7, `30c8d33`): 30/30 cycles nucleate, 2325-2559/2581 nodes rupture per cycle, peak shear 151 MPa, peak normal 105 MPa (compressive), recurrence 1-23 yr | per release | 2026-10-07 | `grep -n slipSense compset/xianshuihe.gmsh.lite/user_defined_params.py` (→ `par.slipSense = -1`) | lars-eriksson |
| Host-side branch/tag protection on `main` — currently unprotected; PR-only + green-required-CI merge policy is now a standing rule (R38) but nothing on the host enforces it yet | P1 | pending owner decision | every visit | 2026-10-07 | (GitHub → Settings → Branches, for `main`) | owner (see PROJECT_RULES.md R38) |
| `subei.gmsh.lite` mixes strike-slip and thrust faults under one `slipSense` value per system — flag only, do not fix | P2 | documented (PR #3, `5f6356e`) | per release | 2026-10-07 | `grep -n "slipSense \* nsmpgp(6" src/interstress.f90` (→ lines 112-113, cited in README) | sophia-okafor (document the limitation) |
| Land the paused 2026-10-04 seed (`git stash@{0}`): `install.sh`/`scripts/case.setup`/`scripts/make_results_bundle.sh` unmarked `&#124;&#124; true` fixes, `test_system/test_conventions.py` checks for R35-R38, `verify_stf.py`/`verify_xianshuihe.py` refactor. The R27a/R35-R38 rule text itself has already been folded into `PROJECT_RULES.md` this pass; only the code/test side of the stash remains | P2 | open | once | 2026-10-07 | `git stash show -p stash@{0}` | lars-eriksson (code), iris-vermeulen (test coverage) |
| Rename `compset/` / `test_system/` to EQdyna's own `case_input/` / `testsys/` naming | P2 | deferred by owner | — | 2026-10-07 | n/a | owner |
| GitHub Support purge request for commit `46a55a7` (private fault-system experiment; already force-removed from this repo's history 2026-10-06) | P3 | pending owner decision | — | 2026-10-07 | n/a | owner |
| Delete the superseded 2026-10-01 results bundle under `work/bundles/saf_modelA_stf_20261001{,.zip}` (16 GB zip + folder, superseded by the `20261005` bundle) | P3 | pending owner decision | — | 2026-10-07 | `du -sh work/bundles/saf_modelA_stf_20261001*` (→ 22G + 16G) | owner |
