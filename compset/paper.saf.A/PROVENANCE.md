# Provenance of the published Liu et al. (2022) artifacts

Paper: *Observation-Constrained Multicycle Dynamic Models of the Southern
San Andreas and the Northern San Jacinto Faults* — DOI
[10.1029/2021JB023420](https://doi.org/10.1029/2021JB023420).

**What this records.** The loading file that produced the published Model A
results is not cleanly present in either public archive. This is the evidence
for that claim, and it is why `compset/paper.saf.A/README.md` carries a
provenance caveat on `Rate_direction.txt`.

**Paths below.** `zenodo_software/` and `pangaea_results/` are fetched by
`bash scripts/fetch_published_reference.sh` (DOIs, md5-pinned) — run that first
if you want to follow the trace yourself. `source_local_recovered/` is a
paper-era tree found on disk by D. Liu and **was never published**, so it is
not fetchable; where the argument below depends on it, that is stated.

The paper published two archives:

| Subdirectory | Source | DOI | Contents |
|---|---|---|---|
| `zenodo_software/` | Zenodo (record 5823021 = GitHub tag v1.0.0) | `10.5281/zenodo.5823021` | Executable `eqdyna2d_2.0.2`, a minimal input set, MATLAB post-processing scripts, raw CFM5.2 fault traces + MATLAB geometry-prep scripts. **No Fortran source.** |
| `pangaea_results/` | Pangaea (dataset 940262) | `10.1594/PANGAEA.940262` | Archived Model A/B/C simulation outputs plus a `mesh/` subdir of template/preprocessing files. **The `mesh/` inputs are NOT what the compiled binary consumes at runtime** — see "Rate_direction.txt trace" below. |
| PDF | Wiley AGU Publications | `10.1029/2021JB023420` | The paper. |
| `source_local_recovered/` | **Not published** — recovered locally | — | Fortran source tree for `eqdyna2d_2.0.2`. Found on disk by D. Liu after the fact; colocated here for development reference. Binary md5 matches Zenodo's. See its `PROVENANCE.md`. |

## Model A / B / C parameters (Table 1 in paper)

| Model | `fric_fs` | `fric_fd` | `fric_fv` (= f_r) | Effective η_min (Pa·s) |
|---|---|---|---|---|
| A (preferred) | 0.5 | 0.465 | 0.49 | 7.0×10²¹ |
| B | 0.3 | 0.265 | 0.29 | 4.2×10²¹ |
| C | 0.7 | 0.665 | 0.69 | 12.0×10²¹ |

`fric_fini` is **not listed in the paper**. Archived inputs show fini for A (0.4) and B (0.25); **fini for Model C is not archived** and must be recovered from the archived initial-stress output in `pangaea_results/work_vis12_fs0.7/totalop.txt1` (first-cycle rows give `σ_shear/|σ_normal| ≈ fric_fini`).

The per-model "minimum viscosity" (7.0 / 4.2 / 12.0 ×10²¹) in the paper's Table 1 is not a direct input. FE_Global.txt sets `ant0 = 8.4×10²¹ Pa·s` (shared across models); at runtime `interstress.f90` uses `ant = ant0` directly. The per-model η values are an effective proxy following the paper's `η_min = -fs·σₙ/γ_max` relation; changing `fric_fs` alone produces the paper's per-model behavior.

## Rate_direction.txt — runtime trace and unit convention

This file is referenced in four `.f90` sources, but only two get compiled into `eqdyna2d_2.0.2` (per `source_local_recovered/code/makefile` OBJ list):

**Runtime consumers (compiled):**
- `eqdyna2d.f90:65–66` opens and reads into global `rd(2, maxftnode*ntotft)`; no post-read scaling.
- `interstress.f90` is the **only place** `rd` is used during multicycle runs. Formula (`line 57–59`):
  ```
  ant = ant0                            ! from FE_Global.txt (= 8.4e21 Pa·s)
  rs  =  rd(1) * cos(2θ) * ant          ! asymptotic shear stress
  rn  = -rd(1) * sin(2θ) * ant          ! asymptotic normal stress
  ```
  `rs`/`rn` are stress values (Pa). Dimensional analysis: `rs = rd(1) · ant` with `ant` in Pa·s requires **`rd(1)` in s⁻¹**. For physical ~50 MPa stress, `rd(1)` ≈ 1e-15 to 1e-14 s⁻¹.

**Standalone programs NOT compiled into the binary** (both are `program ...`, not subroutines, and not in the makefile):
- `strbld.f90`, `strini.f90` — use a different formula, `rs = (str/450)·rd(1)·cos(2θ)·ant`, with `rd(1)` expected in **nrad/yr** (literal values ~100–200). These exist as preprocessing tools; the main cycle driver never calls them.

**Implication for the three Rate_direction.txt variants in this tree:**

| File | Col 1 values | Consistent with `interstress.f90` at runtime? |
|---|---|---|
| `source_local_recovered/input/Rate_direction.txt` | ~6e-15 (rad/s, varying) | ✓ yields ~54 MPa; this is what the binary needs |
| `zenodo_software/dunyuliu/.../input/Rate_direction.txt` | uniform `160 13` (nrad/yr literal) | ✗ yields ~1.3e24 Pa — unphysical |
| `pangaea_results/mesh/Rate_direction.txt` | varying 0–200+ (nrad/yr literal) | ✗ yields ~10²⁴ Pa — unphysical |

**Conclusion:** The "160 13" / "100–200" literal format in the Zenodo and Pangaea `mesh/` archives matches the **unused** `strbld`/`strini` formula, not the binary that actually ran. Those files appear to be legacy or preprocessing artifacts that were archived but never fed to the compiled binary. The only Rate_direction.txt we have in a format consistent with the runtime consumer is the one in `source_local_recovered/input/` (the locally-recovered tree). This is uncomfortable — the authoritative paper-run loading file is therefore **not cleanly in either public archive**; the closest match is the recovered-locally file.

## Reproduction recipe

1. Start from `source_local_recovered/input/` for the seven input files (`FE_Fault_Geometry.txt`, `FE_Global.txt`, `FE_Model_Geometry.txt`, `Rate_direction.txt`, `x{1,2,3}_1.txt`). `Rate_direction.txt` here is in the rad/s convention that the compiled binary consumes.
2. Copy those files plus `zenodo_software/dunyuliu/.../exe/eqdyna2d_2.0.2` into a run directory.
3. Edit only `FE_Global.txt` lines 9–12 per the Model A / B / C table to select the model; line 12 (`fric_fini`) for Model C is unknown and must be recovered from archived `totalop.txt1`.
4. Set `icstart, icend` on line 7 to the paper's cycle counts (A: 4000, B: 4000, C: 4680) or to segments matching the archived `totalop.txt{icstart}` boundaries (A: 1, 1256, 2929; B: 1, 154, 1378, 2888; C: 1, 1436, 3284).
5. Run the binary directly: `./eqdyna2d_2.0.2`. No `case.setup`, no `meshgen.py` — those belong to the post-2.0.2 workflow.

Ground-truth output to diff against lives in `pangaea_results/work_vis7_fs0.5/` (A), `work_vis4.2_fs0.3/` (B), `work_vis12_fs0.7/` (C).

## Inconsistencies identified

### Within the Zenodo archive
- **README ↔ directory layout.** Zenodo README says "With input files located in the same directory [as `./exe/`], `eqdyna2d_2.0.2` should run." Input files actually live under `./input/`; the binary needs files from both `exe/` and `input/` colocated to run.
- **README claims A/B/C reproducibility but only ships B.** `input/FE_Global.txt` carries Model B parameters. Running A or C requires hand-editing friction parameters; the README does not mention this or give the Table 1 values.
- **`input/Rate_direction.txt` cannot be what produced the paper.** Its literal `160 13` values would yield ~10²⁴ Pa stress under the compiled binary's `interstress.f90` formula. The file matches the unused `strbld`/`strini` formula instead.

### Within the Pangaea archive
- **`mesh/Rate_direction.txt` cannot be what produced the archived outputs.** Same issue: literal values 0–200+ under `rs = rd(1)·ant` with `ant=8.4e21` give unphysical stress. The paper's actual per-run loading file is not cleanly archived in Pangaea.

### Between Zenodo `input/` and Pangaea `mesh/`
| File | Zenodo `input/` | Pangaea `mesh/` | Verdict |
|---|---|---|---|
| `FE_Global.txt` | Model B params | Model B params | content identical |
| `FE_Fault_Geometry.txt` | identical | identical | ✓ |
| `FE_Model_Geometry.txt` | LF line endings | CRLF line endings | cosmetic only |
| `Rate_direction.txt` | uniform `160 13` | varying col1, col2=13 | both in the wrong format for the runtime binary (see trace above) |
| `x{1,2,3}_1.txt` | identical | identical | ✓ |
