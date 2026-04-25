# Changelog

All notable changes to EQdyna.2Dcycle. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows [SemVer](https://semver.org/).

## [2.0.7-rc6] - 2026-04-25

Pre-release of 2.0.7 (iteration on rc5). Drops the legacy test.*
compsets, fixes the paper.saf.A binary pin, and adds a paper.saf.A
README so all four compsets are documented.

### Removed
- `compset/test.saf`, `compset/test.subei`, `compset/test.gulang`:
  legacy gmsh test cases superseded by the `.gmsh.lite` versions.
- `compset/paper.saf.A/buildLoadingInput.py`: orphan after test.saf
  removal (only consumer was test.saf input generation).
- `test_system/reference.results/test.subei/`: no longer applicable;
  verify.test.py warns (does not fail) on missing reference, so the
  reference for `subei.gmsh.lite` can be regenerated at leisure.

### Fixed
- `compset/paper.saf.A/user_defined_params.py`: dropped
  `par.exe = 'run_eqdyna2d_2.0.3'` override that pinned the demo to a
  binary no longer in `bin/`. Now falls through to
  `defaultParameters.py` which auto-tracks the repo `VERSION` file,
  matching the other three compsets.

### Added
- `compset/paper.saf.A/README.md`: documents mesh mode (C_mesh=2),
  loading source (`Rate_direction.txt`), run instructions, and
  cross-link to `saf.gmsh.lite` as the C_mesh=3 companion.

### Changed
- `test_system/smoke.py`, `test_system/testNameList.py`: switched
  smoke compset from `test.subei` to `subei.gmsh.lite`.
- `test_system/test_all_quick.py`: dropped the three `test.*` entries;
  COMPSETS now lists only `saf.gmsh.lite`, `subei.gmsh.lite`,
  `gulang.gmsh.lite`. All 3 pass at icend=10 in ~3.5 min.
- `CLAUDE.md`: dropped legacy test.* mentions from the compset
  inventory and reference-cases list.
- `example_workflow.sh`: added a comment block listing the four
  compsets and explaining the optional `python3 meshgen.py` step
  needed for C_mesh=3 cases.
- `scripts/case.setup`: generated `run.sh` now invokes the binary as
  `./<exe>` (case-dir-local) instead of bare `<exe>`. Removes the
  hidden dependency on `install.sh` having been sourced in the caller
  shell to put `bin/` on PATH.

## [2.0.7-rc5] - 2026-04-25

Pre-release of 2.0.7 (iteration on rc4). Adds the `gulang.gmsh.lite`
compset, generalises the post-processing scripts to handle arbitrary
`ntotft`, and ships a quick smoke-test runner.

### Added
- `compset/gulang.gmsh.lite/`: gmsh-meshed (C_mesh=3) Gulang fault
  system (ft1..ft5, all left-strike vertical, uniform-x loading).
  ftVis = 5e21 Pa·s (1e21 too low for nucleation, 8.4e21 gives tensile
  normal at large strikes). ft3 trimmed to drop a sliver-element
  overlap with ft4. First M=7.06 system-spanning event reproduced
  cleanly; ~58 events in 250 yr at icend=4000.
- `compset/gulang.gmsh.lite/preprocess_control_points.py`: resamples
  the dense raw `.gmt` control-point digitisation (412 pts, median
  ~30 m spacing — well below mesh dx=0.3 km) down to a SAF-like
  39-pt set at >=3 km arc-length spacing via per-fault smoothing
  splines (UnivariateSpline on x(s), y(s)). Originals preserved as
  `*.gmt.txt.orig`; idempotent.
- `scripts/meshGenLib.py decimateControlPoints(xs, ys, dx)`: thin a
  control-point series so consecutive kept points are >=dx apart in
  arc length. Avoids spline derivatives tracking sub-mesh jaggedness.
- `scripts/meshGenLib.py plotFaults(...)`: map-view fault layout
  sanity-check plot, colour-coded by ftType with optional loading
  arrow. Saved by saf/subei/gulang.gmsh.lite meshgen.py to
  `aPlots/faults.png`.
- `test_system/test_all_quick.py`: quick smoke-test runner. Builds,
  then for every compset *except* paper.saf.A: runs at icend=10,
  waits for completion, runs the full post-processing suite (catalog
  + b-value + rupture-dynamics + Figure 4).

### Changed
- `scripts/saf_result_utils.py load_fault_blocks`: auto-detects ntotft
  from `meshGeneralInfo.txt` line 2 and falls back to generic names
  `ft1..ftN` when the count differs from the SAF triple. Previously
  hardcoded to 3 SAF blocks, silently mis-slicing totalop rows for
  non-SAF cases.
- `scripts/plot_event_slips_overtime_fig4.py`: per-fault colours now
  auto-fill from `tab10` for non-SAF blocks; SAF colour mapping
  preserved when names match.
- `scripts/meshGenLib.py uniformXLoadingAngle`: sign convention
  corrected to `phi = alpha = atan2(ty, tx)` (paper convention
  `phi = fault_strike - max_shear_direction`, with max-shear along +x
  so max_shear_direction = 0). Was `-atan2(ty, tx)`.
- `scripts/meshGenLib.py writeFilesForEQdyna`: `nsmp/nsmpTanLen/
  nsmpGeoPhys` row counts generalised from hardcoded `x3` to `xnFt`,
  enabling >3 fault systems (gulang has 5).
- `scripts/case.setup`: generated `run.sh` now runs
  `plotRuptureDynamics` BEFORE moving outputs to `aRawSimuData/`
  (was after; the move broke per-cycle plotting).
- `compset/saf.gmsh.lite/meshgen.py`,
  `compset/subei.gmsh.lite/meshgen.py`: call `plotFaults(...)` after
  geometry setup.
- `README.md`: added Outputs and Post-processing sections (file table,
  script table, catalog/analyze quick-start).
- `test_system/smoke.py`: poll `run_*.log` for "Job finished" with a
  600s timeout instead of assuming `bash run.sh` is synchronous (it
  has been backgrounded via nohup since rc1, which silently broke
  the gating). Also switched the build flag from `-m macos` to
  `-m ubuntu` so the gate works on Linux release hosts.

## [2.0.7-rc4] - 2026-04-24

Pre-release of 2.0.7 (iteration on rc3). Adds `subei.gmsh.lite` compset
mirroring the saf.gmsh.lite framework, and the default uniform-x loading
mode for compsets without a paper Rate_direction.txt source.

### Added
- `compset/subei.gmsh.lite/`: gmsh-meshed (C_mesh=3) Subei fault system
  (atf, dxs, sbt) with multi-surface decomposition for connected
  junctions, analytical-spline tangent at every gmsh node, and per-node
  uniform-x loading angle. Bench: 221 events at this writing, M5.0–7.4,
  b≈1.09 (LSQ on [5.0,7.0]), cycle 1 = M7.4 full-system rupture.
- `scripts/meshGenLib.py uniformXLoadingAngle(tx, ty)`: default loading
  model — compression along +x, ftLoadAngle = -atan2(ty,tx) per node.
  Used when no paper Rate_direction.txt is available.
- `scripts/plot_event_slips_overtime_fig4.py`: renamed from
  `plot_saf_figure4_7_8.py` (script content unchanged in this rc;
  generalisation deferred).

### Changed
- CLAUDE.md, scripts/monitor_runs.sh: updated for the rename.

## [2.0.7-rc3] - 2026-04-24

Pre-release of 2.0.7 (iteration on rc2). Adds catalog mode for
plotRuptureDynamics and a catalog-analysis script (b-value, MFD,
nucleation patterns). Strengthens the release.sh docs-lockstep gate.

### Added
- `CATALOG=1 plotRuptureDynamics` mode: fast (no figures), writes
  `aPlots/catalog.csv` with eqId, magnitude, moment, nucleation x/y,
  fault index of nucleation, rupture duration, peak slip per cycle.
  ~1.6s for 438 cycles vs minutes for full plotting.
- `scripts/analyze_catalog.py`: b-value (LSQ on user-windowed
  magnitude range), magnitude-frequency distribution, magnitude vs
  cycle, nucleation along strike vs cycle. Saves
  `aPlots/catalog_analysis.png`.

### Changed
- `release.sh`: docs-lockstep is now the **first hard-fail gate** (was
  last, soft warning). Detects new doc files via `git log --name-only`,
  not just modifications. `--skip-doc-check` to bypass.
- CLAUDE.md compset list now mentions `saf.gmsh.lite` (the lockstep
  doc update missed in rc2).

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
