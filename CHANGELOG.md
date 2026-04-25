# Changelog

All notable changes to EQdyna.2Dcycle. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows [SemVer](https://semver.org/).

## [2.0.7-rc2] - 2026-04-24

Pre-release of 2.0.7 (iteration on rc1). Adds the gmsh-meshed
`saf.gmsh.lite` compset that reproduces paper.saf.A loading on a coarser
mesh, monitor + plot improvements, and self-tracking version banner.

### Added
- `compset/saf.gmsh.lite/` — coarser-resolution (dxy=400m) gmsh compset
  that reproduces Liu et al. 2022 Model A. Self-contained: ships with
  paper.saf.A's `Rate_direction.txt + x{1,2,3}_1.txt`. Verified: cycle 1
  first interval 473 yr (paper: 471 yr), M5.7 NW patch matches paper Fig 4.
- `scripts/meshGenLib.py` helpers: `splineFaultFromControl` (natural cubic
  spline y(x) sampled at dxy, mirrors meshgen1.f90), `loadPaperRateDirAlongArcLen`
  (rebuilds paper per-fault arc-length and γ, φ from x*_1.txt + Rate_direction.txt),
  `plotLoadingInputs` (paper-Fig-2-style γ/φ/wt/vis plot).
- `scripts/monitor_runs.sh` — periodically tails progress + re-plots
  Figure 4 and rupture dynamics for all `paper.saf.A.*` and saflite cases.
- `scripts/plotRuptureDynamics`: `MIN_PLOT_MAGNITUDE` env var (default 6.5)
  with `>=` comparison; skip-existing-plot optimization (set `FORCE_REPLOT=1`
  to override); robust loader for Fortran 3-digit-exponent output.
- `src/makefile` cpp macro `EQDYNA_VERSION` injects the VERSION file
  string into the runtime banner — no more hardcoded `2.0.0`.

### Changed
- `src/interstress.f90` C_mesh=3 branch: comments document the
  paper-faithful encoding (ftVis = ant0·str/γ per node) where
  `ant = ftVis·450/ftLoadWt = ant0·str/γ` matches the C_mesh=2 paper formula.
- `compset/saf.gmsh.lite/meshgen.py`:
  - Single SAF control file (no more ssaf1+ssaf2+bridge split).
  - Per-node ftLoadAngle, ftVis interpolated from paper Rate_direction.txt
    instead of uniform constants.
  - C¹-smooth tangent at every gmsh fault node via analytical spline
    derivative (matches meshgen1.f90 lines 437-442).
- `scripts/case.setup`: C_mesh=3 branch now exports
  `GFORTRAN_UNBUFFERED_ALL=1` and `OMP_NUM_THREADS=${...:-1}` (was
  C_mesh=2 only); 1-thread default replaces previous 16-thread default.
- `scripts/defaultParameters.py exe`: walks up directory tree to find
  `VERSION` file (handles case-dir invocation where `__file__/..` doesn't
  resolve to repo root).
- `install.sh`: `mv` (not `cp`) binary into `bin/` to leave `src/`
  source-only; `rm src/*.o src/*.mod` after build to keep `src/` tidy.
- Welcome banner: cleaned attribution (UT Austin + Texas A&M),
  contact email updated, mentions `nsmpGeoPhys.txt` for C_mesh=3,
  `cylce`→`cycle` typo fixed.

### Fixed
- `cp src/run_eqdyna2d_* bin/` failed with "Text file busy" when a
  running simulation had the binary mmap'd → install.sh now `rm -f` the
  old binary first (unlinks the inode; running procs keep their copy).

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
