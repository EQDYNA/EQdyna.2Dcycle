# Changelog

All notable changes to EQdyna.2Dcycle. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows [SemVer](https://semver.org/).

## [2.0.7-rc1] - 2026-04-24

Pre-release of 2.0.7. First release on top of published `v2.0.6`. Surface-level changes (install, README, compsets, OMP, release process) are stabilized; soliciting feedback before final 2.0.7 tag. Test-system reference results not yet refreshed against the new compsets.

### Added
- `paper.saf.A` compset reproducing Liu et al. (2022) JGR Solid Earth Model A (C_mesh=2, smooth per-node `Rate_direction.txt` recovered from paper-era local archive).
- `example_workflow.sh` — one-step demo (install + create case + setup + launch).
- `VERSION` file at repo root drives binary name (`bin/run_eqdyna2d_<VERSION>`).
- `test_system/smoke.py` — quick build + 1-cycle sanity check.
- `test_system/README.md` — usage and "add a new test case" guide.
- `CHANGELOG.md`, `release.sh` — single-source release process.
- OMP parallelization with `$OMP ATOMIC` in `qdct3` and `hrglss`; per-op timing (`t_vd`, `t_qdct3`, `t_hrglss`, `t_faulting`, `t_brhs`) printed every 300 steps; date/time stamps in driver log. ~5.7× speedup at 8 threads on paper.saf.A.
- Rate-dependent viscosity in `interstress.f90` C_mesh=2 branch: `ant = ant0 * str / rd(1)` (paper formula, replaces uniform fallback).
- `create.newcase` flags: `--work_dir`, `--compset`, `--force`, `--list`.
- Auto-detect OS in `install.sh`; installs python deps on both Linux and macOS.
- `run.sh` exports `GFORTRAN_UNBUFFERED_ALL=1` and `OMP_NUM_THREADS=${OMP_NUM_THREADS:-16}`.

### Changed
- Renamed `install.eqdyna.2dcycle.sh` → `install.sh`.
- Renamed `test/` → `test_system/` (avoids Python stdlib collision); `testAll.py` → `test_all.py`.
- `test_system/verify.test.py` now exits non-zero on any failure; `test_all.py` propagates verify status.
- `install.sh` rebuilds binary safely when an old copy is in use (rm-then-cp pattern).
- README rewritten: install autodetect form first; Run section shows both `example_workflow.sh` and the manual 4-command flow.

### Removed
- `makeExecutables.sh` (one-line chmod helper, redundant).
- `test.sh` (overlapped with `example_workflow.sh` and `test_system/`; name clashed with `test_system/`).
- Per-release `release_notes_*.md` files (replaced by this single `CHANGELOG.md`).
- `docs/` directory (moved to untracked `archive/`; release notes consolidated here).

### Fixed
- Velocity/displacement update, `qdct3`, `hrglss`, and `brhs *= alhs_inv` are thread-safe under OMP (verified bit-identical results across 1/2/4/8 threads).
- `install.sh` no longer fails with "Text file busy" when a running simulation has the binary mmap'd.
- Driver per-step print format string (16 descriptors / 17 args bug).

---

History prior to v2.0.3 was tracked in per-release `release_notes_*.md` files; those are no longer maintained. See `git log --all --tags` for the historical record.
