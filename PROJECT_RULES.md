# Project rules — EQdyna.2Dcycle

Rules earned from actual incidents, not style preferences. Each states the
convention, why it exists, and how it is enforced. Every rule with an `R#`
id has a check in `test_system/test_conventions.py`; run it with
`python3 -m test_system.test_conventions` — fast (no binary, no meshing),
meant to run on every change.

---

## Mesh file conventions

### R1 — Never subtract 1 from `fac.txt` / `nsmp.txt` unconditionally
`meshgen.py` writes these **0-based**; some older cases are 1-based. Shift
only when the minimum index is 1:
```python
f = np.loadtxt('fac.txt', dtype=int)
if f.min() == 1:
    f = f - 1
```
An unconditional `- 1` silently scrambles every cell — it once inflated a
computed mesh area to 19.3M km² against a 139k km² domain and drew a mesh
of random crossing lines, with no error raised. `plotMeshNearFault.py` and
`checkMeshQuality.py` already guard this; copy from them.

### R2 — `nsmp.txt` is padded per fault; slice by `nfn[i]`, don't filter on `id > 0`
Each fault's block is `maxftnode` long but only the first `nfn[i]` rows are
its nodes: `ids = nsmp[i*maxftnode : i*maxftnode + nfn[i], 0]`. Filtering
`ids > 0` drops node 0 (real, 0-based), and subtracting 1 from a padding `0`
wraps to the last vertex, drawing a spurious line across the domain.

### R3 — Fault polylines must be strictly increasing in x
`meshgen.py` fits a `CubicSpline(x, y)`, which requires strictly increasing
x; merging polylines that share an endpoint leaves a duplicated x and the
spline rejects it — `seg0`/`seg4` of the Xianshuihe trace share an endpoint
at exactly 0.000 km. Run `scripts/check_fault_geometry.py` before meshing
any new fault system.

### R18 — Classify fault sides by the fault normal, never by global dy
`judgeElemDirect` used to decide above/below from the **global sign of dy**
between on-fault and off-fault centroids — valid only on gently-dipping
faults. On a steep section both sides have dy > 0, both get labelled
"above", and `replaceMasterWithSlaveNodes` swaps all of them to the slave,
leaving the master node with no cell. Use the local fault tangent and the
sign of the cross product instead. This orphaned one master node on steep
Xianshuihe ft4, and silently misclassified 3 cells (0.022%) in
`subei.gmsh.lite` — never orphaned, never noticed, traction applied through
the wrong side. An orphaned split node has no cell to transmit traction;
`checkMeshQuality` reports orphans, run it on every new mesh.

### R20 — The mesh export must not drop elements
`meshGenLib` writes only quads to `fac.txt`; where gmsh cannot recombine a
region it leaves triangles, discarded without a word, leaving holes in the
mesh. The Xianshuihe mesh had 34602 quads + 10 triangles in the `.msh` but
34602 rows in `fac.txt` — the lost triangles sat on the two orphaned split
nodes and were badly shaped (min angle 15°). `checkMeshQuality` couldn't
see it, since it reads `fac.txt`, and reported `angle<20 deg = 0` after
those elements were already dropped. Compare `fac.txt`'s element count
against the `.msh`'s 2D element count on every mesh; see issue #1.

### R21 — Every mesh failure becomes a mechanical guard
Fix the mesh, then add a check that makes the same defect impossible to
ship unnoticed — a **hard failure**, non-zero exit, not a warning.
`checkMeshQuality.py` hard-fails on:

| check | why it exists |
|---|---|
| orphaned master/slave split node | no element to transmit traction |
| element count vs the `.msh` | `fac.txt` writes quads only; dropped triangles leave holes |
| triangles present at all | gmsh could not recombine, usually a sub-element fault gap |
| interior angle < 20 deg | sliver |
| interior angle > 160 deg | near-flat corner; the report used to check only the minimum, so 164.8 deg passed |
| aspect ratio > 10 | stretched cell |
| degenerate cell | zero area |
| mixed master/slave cells | split-node bug |

Every one came from a defect that reached a mesh and was missed (see the
orphan hunt above). Thresholds are named constants at the top of
`checkMeshQuality.py` (`MIN_ANGLE_DEG`, `MAX_ANGLE_DEG`, `MAX_ASPECT`,
`MIN_FAULT_GAP_ELEMENTS`).

### R28 — Fortran must never enumerate faults by hand; loop over `ntotft`
Fixed hand-written branches (if/elseif chains, or an index summing literal
`nfnode(1) + nfnode(2) + nfnode(3)` terms) break the moment `ntotft` exceeds
that number, silently if there is no `else`. `interstress.f90` excluded
fault tips and searched for nucleation on faults 1-3 only, leaving
`ntotft > 3` unable to host an earthquake there. `faulting.f90` resolved
the nucleation point with four if/elseif branches and no else; a nucleation
on fault 5+ left `ift0`, `xcoor0`, `ycoor0` uninitialised, the
forced-rupture patch never applied, and every interseismic period collapsed
to the 1-year floor — endless "cycles" with no earthquakes and nothing in
the log to say why. Both were silent on 3-fault cases, surfacing only on
the 5-fault (`gulang`) and 7-fault (`xianshuihe`) compsets. Mechanical
check: `test_no_hardcoded_fault_counts_in_fortran` bans three or more
summed literal `nfnode(<digit>)` terms in `src/*.f90` — it does not catch
an unbounded if/elseif chain with no arithmetic, a Tier-3 gap.

---

## New fault system: decisions that cannot be inferred

### R4 — Run `check_fault_geometry.py` before wiring a new compset
It reports the recurring interventions: faults too short to resolve
(`RESOLUTION`), faults sharing a node (`SHARED NODE`), faults closer than
one element (`SUB-ELEMENT`), duplicate x (`DUPLICATE`) — nothing is
auto-fixed, each is a modelling decision. Both Xianshuihe problems (a 15 km
fault, two faults meeting at 0.000 km) were found only because distances
happened to be checked by hand.

### R5 — Don't merge polylines to tidy geometry; merge only with a stated reason
Merging erases structure — step-over density is a physical quantity (Qiao
et al. 2022 use it as a maturity proxy), so concatenating polylines
destroys the thing being modelled. Legitimate reasons, both used for
Xianshuihe: the junction is far below the step-over cut-off (a digitisation
break, not a structure), or a fault is too short to resolve. Record the
reason in the compset README.

### R6 — Parallel strands cannot be merged into a chain
End-to-end joins collapse into one fault; side-by-side strands need
multi-surface topology (`subei.gmsh.lite`), not chain topology
(`gulang.gmsh.lite`).

### R19 — Every compset parameter must be assigned to `par`
`case.setup` reads attributes off `par`, so a bare `name = value` in
`user_defined_params.py` sets a module-level local nothing ever reads.
`fric_fini = 0.45` was missing the prefix in `gulang`, `subei` and
`xianshuihe`; it never had any effect, unnoticed only because the value
happened to equal the default in `defaultParameters.py` — change the
number and nothing would happen, silently.

---

## Loading conventions

### R7 — `ftLoadMaxShear` is a strain rate that sets a stress, not a slip rate
`interstress.f90` computes `rs = rd·cos(2θ)·ant`, an asymptotic **stress**.
The long-term slip rate is emergent, so a published slip rate cannot be
substituted directly, and inverting `γ̇ = V/2W` against the mesh width gives
nonsense (~570 mm/yr for the SAF). Calibrate: scale from the SAF anchor
(`1.427e-14 s⁻¹` → ~38 mm/yr peak), run, measure with `plot_saf_figure3.py`,
rescale linearly.

### R8 — `ftLoadAngle = -999` is a placeholder, not a measurement
It makes the code assume loading is along x and derive θ from the fault
tangent — fine where the rotated frame aligns with plate motion, elsewhere
it discards the real orientation of the strain field. A proper strain-rate
map supplies both the max shear magnitude and its direction.

### R9 — Domain extension and `term` are coupled
`ext = 10 + totalSimuTime*vp/1e3` km keeps boundary reflections from
reaching the fault before the rupture ends. Shrinking `ext` for a cheaper
mesh breaks that guarantee — at `ext = 100 km` a P-wave round trip is ~33 s
against a 200 s run. Acceptable for meshing and geometry checks; state it
in the file and restore before trusting rupture dynamics.

---

## Reproducing published results

### R10 — Parity is against the reference implementation running, not its printed output
Gate a port on executing the original (MATLAB is on PATH) against the
published data — printed tables are rounded and can be wrong. Published
Table 2 prints 49 events at Frazier Mountain where its own script on its
own published data computes 52; matching the paper would have meant
breaking parity with the code that made it.

### R11 — Reproduce reference bugs deliberately, and say so
MATLAB's `load()` reads gfortran's dropped-E `0.8384675-101` as `0.8384675`
— truncation, not scientific notation. Reproducing published numbers
requires truncating too. Document any such choice at the point of use;
never silently "fix" a reference implementation while claiming parity.

### R12 — State the magnitude convention wherever magnitudes are compared
`plotRuptureDynamics` uses the model's own `mu = rou*vs^2 = 3.204e10 Pa`;
the paper's Figures 6 and 9 use `3500^2*3000 = 3.675e10 Pa`. Every catalog
magnitude and b-value therefore sits 0.04 Mw below the paper's scale.

### R13 — Constants derived from an external chain must record their derivation
`OBSERVED_EQDYNA_X_KM` was wrong by up to 147 km, misassigning a third of
the SSAF nodes to the wrong observed-rate band in every Figure 3 ever
produced, with no comment saying where it came from — nothing flagged it.
Derive, verify against the source chain, and record how.

### R14 — Hand-picked cycle ids do not transfer between runs
Sequences are chaotic; selecting events by cycle id from another run
mislabels unrelated events — Figure 9 stamped "1857-like" on an M5.7 that
broke 8 km of fault. Select by rupture footprint, label only on a match.

---

### R30 — Results reproduce only at a fixed `OMP_NUM_THREADS`
The OpenMP reductions in `qdct3` and `hrglss` accumulate into `brhs` under
`$OMP ATOMIC`, so the summation order depends on thread scheduling and the
result differs in the last bit between thread counts. Earthquake sequences
are chaotic: that difference amplifies to a visibly different catalogue
within about four cycles. Any run being compared against a stored reference
— a benchmark, a regression check, a published result — must use the same
thread count as the reference, and the thread count belongs in the record.

Verified on v2.1.0 against both meshing paths: at `OMP_NUM_THREADS=1` both
`paper.saf.A` (C_mesh=2) and `saf.gmsh.lite` (C_mesh=3) reproduce their
stored `totalop.txt1` with max absolute difference **0.0** over 4 and 10
cycles; at 2 threads the same binary matches for 3-4 cycles then diverges
(cycle 5 interval 53 yr against the reference 128 yr). Do not read such a
divergence as a regression before ruling this out.

Mechanical check: `test_run_sh_sets_omp_threads` — the generated `run.sh`
must set `OMP_NUM_THREADS` explicitly rather than inheriting whatever the
caller's environment happens to hold.

## Case I/O

### R15 — Read case outputs from the case dir **or** `aRawSimuData/`
`run.sh` moves `totalop.txt*`, `cyclelog.txt*`, `interval.txt*` (and
`binaryop` for C_mesh=2) into `aRawSimuData/` when a run finishes; a script
that only looks in the case dir breaks on every completed run. Search both,
never sum both — double-counting `interval.txt` inflated the simulated
duration.

### R16 — `run.sh` must not delete `binaryop` before launching
It is the restart state; deleting it made the documented
restart-from-binaryop procedure impossible for C_mesh=2 cases.

### R29 — A plot script must never render an empty result without saying so
A figure with axes, labels and a legend but no data reads as a finished result.
`plot_event_slips_overtime_fig4.py` did exactly that twice: once when its
SAF-tuned `--threshold 1.0` m filtered out every event of a sub-metre
catalogue, and once when a requested window started past the last event —
event times begin at the first event, so the last one sits at
`sum(intervals) - intervals[0]`, well short of total simulated time on a run
with a long first interseismic period. Both saved a plausible-looking PNG and
exited 0. A script that draws nothing must say what it filtered and why, and
exit non-zero. Mechanical check: `test_fig4_refuses_empty_output`.

### R17 — Detach the run wrapper
`setsid` where available, so a process-group kill cannot orphan the binary
and skip the post-run steps. A `paper.saf.A` run lost its wrapper
mid-flight: the binary finished all 4000 cycles, but nothing was plotted,
nothing moved to `aRawSimuData/`, and no `Job finished` marker was written.

---

## Release process

`release.sh` was removed on 2026-08-28 — the release process is manual
again. The gates it encoded are restated here so they are not lost with the
script. R24/R26/R27 have mechanical checks in
`test_system/test_conventions.py`; the rest are procedural, run by hand at
release time.

### R22 — Refuse to release from a dirty tree
`git status --porcelain` must be empty before tagging — the point of a tag
is that `git checkout v$VERSION` reproduces exactly what was tested, and
uncommitted changes break that silently. Apply: `git diff --quiet &&
git diff --cached --quiet` before tagging. Procedural — the tree is
legitimately dirty during normal development, so this is a release-time
gate, not a standing check.

### R23 — Docs lockstep at release: README.md / CLAUDE.md / any compset README must have been touched since the previous tag
`git log --name-only <prev-tag>..HEAD -- README.md CLAUDE.md
test_system/README.md work/README.md compset/*/README.md` must be
non-empty. This is starter rule 11 ("docs move with the code") applied at
the release boundary, because that's where it was actually lost —
individual commits can each look fine while the release as a whole drifts
from what its CHANGELOG claims. Apply: run that `git log --name-only`
command against `git describe --tags --abbrev=0` before tagging.
Procedural — depends on the previous tag and commit range, neither fixed
at a given commit.

### R24 — CHANGELOG must have a real body under `## [<VERSION>]`, not just `[Unreleased]`
`CHANGELOG.md` must have a `## [<VERSION>]` heading (VERSION from the root
`VERSION` file) with at least one non-blank body line before the next `##
[` heading or EOF. A release with nothing under its heading means either
the work was never described, or the heading was added without moving the
`[Unreleased]` notes under it — the same failure either way. Mechanical
check: `test_changelog_has_body_for_version`.

### R25 — Run `test_system/smoke.py` before tagging
The last gate before a tag; it must exit 0 — the cheapest thing that would
have caught a release built against a broken build before the tag, not
before someone using the tag, finds out. Apply: `python3
test_system/smoke.py` immediately before `git tag`. Procedural — running a
test suite is not itself a repo-state check `test_conventions.py` can
assert.

### R26 — The release tag message is the CHANGELOG section body — extract it correctly
The naive awk range `/^## \[$VERSION\]/,/^## \[/` is wrong: the start line
also matches the end pattern (both are `## [...]` headings), so the range
closes on the same line it opens and yields nothing. Every tag through
`v2.0.7-rc7` carries an empty annotated-tag message because of exactly this
bug — `release.sh` ran the broken range for its entire life and nobody read
a tag message closely enough to notice. Apply: skip the heading line
itself, then stop at the *next* heading:
```awk
awk -v hdr="## [${VERSION}]" '
    index($0, hdr) == 1 { found = 1; next }
    found && /^## \[/ { exit }
    found { print }
' CHANGELOG.md
```
Reject an empty (whitespace-only) result — R24's check, applied to the
exact text that would ship as the tag message. Mechanical check:
`test_changelog_section_extraction_returns_text` runs this exact
extraction (in Python, same two-condition logic) against the current
`VERSION` and asserts it is non-empty.

### R27 — `VERSION` is the single source of the release version
The root `VERSION` file is a single-line semver (optionally with a
pre-release suffix, e.g. `2.0.7-rc8`), read by `src/makefile` to name the
binary `run_eqdyna2d_$(VERSION)`, and is what `v$(cat VERSION)` tags
against. Nothing else — not a hardcoded string, not a compset README —
states the version independently: two sources diverge the moment one is
bumped and the other isn't; `paper.saf.A` was found hardcoding `par.exe =
run_eqdyna2d_2.0.3` instead of reading `VERSION` (fixed in b3ea71b).
Mechanical check: `test_version_file_is_single_semver_line`.
