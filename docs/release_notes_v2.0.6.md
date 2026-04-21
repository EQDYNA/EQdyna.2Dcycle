# EQdyna.2Dcycle v2.0.6 Release Notes

**Release Date**: April 21, 2026

## Summary of Scope

This release stabilizes the newer `C_mesh = 3` Python/Gmsh workflow, with the main focus on the SAF transition case and the shared case/template machinery that supports `create.newcase -> case.setup -> bash run.sh`.

The release also adds Python plotting utilities for SAF-style figures, fixes `C_mesh = 3` rupture-dynamics magnitude handling, and cleans up release hygiene by excluding generated mesh outputs, caches, and compiled artifacts from version control.

## Files Added, Removed, Renamed, or Cleaned Up

### Added
- `docs/release_notes_v2.0.6.md`
- `scripts/checkMeshQuality.py`
- `scripts/compareMeshGeoPhys.py`
- `scripts/plotMeshNearFault.py`
- `scripts/plot_saf_figure3.py`
- `scripts/plot_saf_figure4_7_8.py`
- `scripts/saf_result_utils.py`
- `src.v2.0.3/` source snapshot and executable tree
- `src/strini.f90`
- `src/strbld.f90`
- `src/TERRArun.slurm`
- `test.sh`
- `compset/test.saf.published/`

### Removed or No Longer Tracked
- `RELEASE_NOTES_v2.0.5.md` moved from repo root to `docs/release_notes_v2.0.5.md`
- `src/catmullrom.f90`
- tracked generated artifacts removed from the index:
  - `compset/test.subei/__pycache__/...`
  - `compset/test.subei/fem_mesh_output/meshWOSplitNode.png`
  - `compset/test.subei/fem_mesh_output/meshWithFaultNodes.png`
  - `compset/test.subei/fem_mesh_output/nsmpGeoPhys.txt`
  - `src/globalvar.mod`

### Cleaned Up
- `.gitignore` now ignores generated case products and build artifacts more aggressively:
  - `compset/**/fem_mesh_output/`
  - `compset/**/__pycache__/`
  - generated `FE_*` case files
  - generated `run.sh`
  - generated mesh/text outputs under `compset/**`
  - `src/*.o`, `src/*.mod`, `src/run_eqdyna2d_*`
  - matching `src.v2.0.3` build artifacts
  - local `.claude/`, `.codex`, and `scripts/__pycache__/`

## Content Updates to Master Documents

- Release-note lineage is now under `docs/`.
- `docs/release_notes_v2.0.5.md` preserves the prior release note after moving it out of the repository root.

## Audit Findings and Fixes

### 1. Generated files were leaking into the repository
Findings:
- tracked Python cache files
- tracked mesh preview images
- tracked generated `nsmpGeoPhys.txt`
- tracked compiled Fortran module output
- untracked local build trees and generated case files

Fixes:
- removed the tracked generated files listed above from the git index
- expanded `.gitignore` to prevent the same classes of files from reappearing in future releases

### 2. `create.newcase` was copying stale generated case files
Findings:
- generated `run.sh`, `FE_Global.txt`, `FE_Model_Geometry.txt`, `FE_Fault_Geometry.txt`, mesh folders, and plot/output directories could bleed from a compset into a fresh case

Fixes:
- `scripts/create.newcase` now skips generated directories and generated case files when creating a new case
- `os.makedirs(..., exist_ok=True)` replaced a narrower directory creation path

### 3. `case.setup` did not cleanly support both legacy and new meshing paths
Findings:
- the generated `run.sh` logic did not clearly distinguish `C_mesh = 2` and `C_mesh = 3`

Fixes:
- `scripts/case.setup` now emits different `run.sh` behavior for:
  - `C_mesh = 2`: legacy internal meshing / `Rate_direction.txt` path
  - `C_mesh = 3`: external Python/Gmsh meshing path
- both paths now generate timestamped log files and launch background jobs consistently

### 4. `C_mesh = 3` SAF tangent/normal handling did not match the intended spline behavior
Findings:
- the SAF development path used a spline-based tangent computation, but it was not aligned with the legacy natural-spline intent

Fixes:
- `scripts/meshGenLib.py`
- `compset/test.saf/meshGenLib.py`

Both now use natural cubic splines for tangent/length evaluation.

### 5. SAF `C_mesh = 3` loading/viscosity mapping was inconsistent with the old variable-`ant` logic
Findings:
- local loading-rate information was being mapped into `ftLoadMaxShear`
- `ftLoadWt` remained constant
- this prevented the intended heterogeneous `ant = ftVis * 450 / ftLoadWt` behavior

Fixes:
- `compset/test.saf/meshgen.py` now maps the interpolated loading magnitude into `ftLoadWt`
- `ftLoadMaxShear` is kept as the constant reference loading rate
- `ftLoadAngle` remains spatially varying
- `ftVis` remains the base viscosity input

### 6. New-case workflow validation was missing for the updated SAF compset
Findings:
- the release needed an end-to-end check that `compset/test.saf` can seed a fresh working case correctly

Fixes:
- validated `create.newcase -> case.setup -> bash run.sh` against a fresh `test.saf` workflow case
- verified that the resulting mesh and first few cycles matched the working `test.saf` reference case

### 7. `plotRuptureDynamics` miscomputed magnitude for `C_mesh = 3`
Findings:
- `nsmpTanLen.txt` lengths for the `C_mesh = 3` path are in km, but the moment calculation was treating them as meters

Fixes:
- corrected moment/magnitude handling in `scripts/plotRuptureDynamics`
- plotting now uses an explicit `Mw > 6.5` save threshold

## Scientific / Workflow Updates Included in This Release

- default pseudo-viscosity in `scripts/defaultParameters.py` updated to `8.4e21`
- `scripts/defaultParameters.py` now defines `exe = 'run_eqdyna2d_2.0.3'`
- `compset/test.gulang` and `compset/test.subei` were updated toward the newer organized mesh-output workflow
- SAF geometry inputs under `compset/test.saf/user_fault_geometry_input/` were updated with the current case content
- several core Fortran files under `src/` were updated as part of the `2.0.3` line, including:
  - `Read_Input_Files.f90`
  - `eqdyna2d.f90`
  - `interstress.f90`
  - `meshgen1.f90`
  - `makefile`

## Remaining Open Issues, Unknowns, or Pending Items

- `C_mesh = 3` generated `run.sh` still calls plain `python3 meshgen.py`; release users need a `python3` environment with `gmsh` installed.
- The remote GitHub repository currently has a wrong version line (`v2.1.0`) that needs to be replaced with `v2.0.6`.
- This release intentionally ignores `archive/`, `docs/`, `misc/`, and `CLAUDE.md` during project audit, except for release-note management under `docs/`.
- There are still many substantive source changes under `src/`; this release note reflects the current final tree but does not attempt to split those changes into smaller historical units.

## Totals or Cost Changes

- No monetary or booking totals apply to this repository.
- No FX assumptions apply.

## Assumptions Used

- Versioning follows the repository’s release-note lineage, not the mistaken remote tag line.
- Therefore the correct next version after `v2.0.5` is `v2.0.6`.
- Release-note management is exempted from the user’s general `docs/` ignore rule because the release workflow explicitly stores release notes there.
