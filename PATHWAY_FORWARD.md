# Pathway forward — EQdyna.2Dcycle

Present-tense status board (invariant 12 / R37): every open issue, to-do and
standing claim, one row, with a priority the work is taken in. State is not
priority — re-prioritise by editing the Priority column only. Blank
"Last checked" = never audited, stays blank. A claim with no command is
remembered, not verified, and is marked so.

| Item | Priority | State | Recheck | Last checked | Command | Route |
|---|---|---|---|---|---|---|
| Release v2.3.0: CHANGELOG `[Unreleased]` holds the `slipSense`/ftType fix; `VERSION` says 2.2.0 but that version was never tagged (latest tag `v2.1.1`) — fold the untagged 2.2.0 content into 2.3.0 rather than tagging 2.2.0 after the fact | P1 | open | per release | 2026-10-07 | `git describe --tags --abbrev=0` (→ `v2.1.1`) vs `cat VERSION` (→ `2.2.0`) | haruto-nakamura |
| `xianshuihe.gmsh.lite` and `gulang.gmsh.lite` are left-lateral fault systems but were simulated right-lateral (solver never read `ftType`, fixed by `par.slipSense`) — set `slipSense = -1` in both, rerun, redo Xianshuihe's target-stress guidance (this session's loading audit: left-lateral resolved shear at 2545/2581 nodes, so the existing 90-100 MPa note is stale) | P1 | open | once, re-verify after rerun | 2026-10-07 | `grep -n slipSense compset/xianshuihe.gmsh.lite/user_defined_params.py compset/gulang.gmsh.lite/user_defined_params.py` (currently no match in either) | lars-eriksson (config + rerun), sophia-okafor (README guidance) |
| Host-side branch/tag protection on `main` — currently unprotected; PR-only + green-required-CI merge policy is now a standing rule (R38) but nothing on the host enforces it yet | P1 | pending owner decision | every visit | 2026-10-07 | (GitHub → Settings → Branches, for `main`) | owner (see PROJECT_RULES.md R38) |
| `subei.gmsh.lite` mixes strike-slip and thrust faults under one `slipSense` value per system — flag only, do not fix | P2 | flagged | per release | 2026-10-07 | `grep -n ftType compset/subei.gmsh.lite/*.py` | sophia-okafor (document the limitation) |
| Land the paused 2026-10-04 seed (`git stash@{0}`): `install.sh`/`scripts/case.setup`/`scripts/make_results_bundle.sh` unmarked `\|\| true` fixes, `test_system/test_conventions.py` checks for R35-R38, `verify_stf.py`/`verify_xianshuihe.py` refactor. The R27a/R35-R38 rule text itself has already been folded into `PROJECT_RULES.md` this pass; only the code/test side of the stash remains | P2 | open | once | 2026-10-07 | `git stash show -p stash@{0}` | lars-eriksson (code), iris-vermeulen (test coverage) |
| Rename `compset/` / `test_system/` to EQdyna's own `case_input/` / `testsys/` naming | P2 | deferred by owner | — | 2026-10-07 | n/a | owner |
| GitHub Support purge request for commit `46a55a7` (private fault-system experiment; already force-removed from this repo's history 2026-10-06) | P3 | pending owner decision | — | 2026-10-07 | n/a | owner |
| Delete the superseded 2026-10-01 results bundle under `work/bundles/saf_modelA_stf_20261001{,.zip}` (16 GB zip + folder, superseded by the `20261005` bundle) | P3 | pending owner decision | — | 2026-10-07 | `du -sh work/bundles/saf_modelA_stf_20261001*` (→ 22G + 16G) | owner |
