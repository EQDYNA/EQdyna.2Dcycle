# Project rules — EQdyna.2Dcycle

Rules earned from actual incidents, not style preferences. Each one states
the convention, why it exists, and how it is enforced. Every rule with an
`R#` id has a check in `test_system/test_conventions.py`; run it with

```bash
python3 -m test_system.test_conventions
```

It is fast (no binary, no meshing) and is meant to run on every change.

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

**Why:** an unconditional `- 1` silently scrambles every cell. It inflated a
computed mesh area to 19.3M km² against a 139k km² domain, and produced a
mesh figure that looked like random crossing lines. It does not raise — it
just makes wrong pictures and wrong areas.

`plotMeshNearFault.py` and `checkMeshQuality.py` already do this; copy from
them, don't reinvent.

### R2 — `nsmp.txt` is padded per fault; slice by `nfn[i]`, don't filter on `id > 0`

Each fault's block is `maxftnode` long but only the first `nfn[i]` rows are
its nodes:

```python
ids = nsmp[i*maxftnode : i*maxftnode + nfn[i], 0]
```

**Why:** filtering `ids > 0` drops node 0, which is a real node under 0-based
indexing. And subtracting 1 from a padding `0` wraps to the last vertex,
drawing a spurious line clear across the domain.

### R3 — Fault polylines must be strictly increasing in x

`meshgen.py` fits a `CubicSpline(x, y)`, which requires strictly increasing
x. Merging two polylines that share an endpoint leaves a duplicated node.

**Why:** `seg0` and `seg4` of the Xianshuihe trace share an endpoint at
exactly 0.000 km. Merging them without dropping the coincident node produces
a repeated x and the spline rejects it.

Run `scripts/check_fault_geometry.py` before meshing any new fault system.

### R18 — Classify fault sides by the fault normal, never by global dy

`judgeElemDirect` originally decided whether a cell sat above or below a
fault from the **global sign of dy** between the on-fault and off-fault
centroids. That is only valid where the fault runs gently. On a steep
section, cells on BOTH sides have dy > 0, so both are labelled "above",
`replaceMasterWithSlaveNodes` swaps every one of them to the slave, and the
master node is left with no cell at all.

Pass the local fault tangent and take the sign of the cross product.

**Why:** this orphaned exactly one master node on the steep NE end of the
Xianshuihe ft4, and silently misclassified 3 cells in `subei.gmsh.lite`
(0.022% of its mesh) that never orphaned and so were never noticed --
their traction was applied through the wrong side of the fault.

An orphaned split node has no cell to transmit traction. `checkMeshQuality`
reports orphans; run it on every new mesh.

### R20 — The mesh export must not drop elements

`meshGenLib` writes only quads to `fac.txt`. Where gmsh cannot recombine a
region it leaves triangles, and those are discarded without a word, leaving
holes in the FE mesh.

**Why:** the Xianshuihe mesh had 34602 quads + 10 triangles in the `.msh`
and 34602 rows in `fac.txt`. The 10 lost triangles sat exactly on the two
orphaned split nodes, and were badly shaped (minimum angles down to 15
degrees). `checkMeshQuality` could not see any of it, because it reads
`fac.txt` -- it reported `angle<20 deg = 0` when five such elements had
already been thrown away.

Compare the element count in `fac.txt` against the 2D element count in the
`.msh` on every mesh. See issue #1.

---

## New fault system: decisions that cannot be inferred

### R4 — Run `check_fault_geometry.py` before wiring a new compset

It reports the recurring interventions: faults too short to resolve
(`RESOLUTION`), faults sharing a node (`SHARED NODE`), faults closer than one
element (`SUB-ELEMENT`), and duplicate x (`DUPLICATE`). Nothing is auto-fixed
— each is a modelling decision.

**Why:** both Xianshuihe problems (a 15 km fault, and two faults meeting at
0.000 km) were found only because distances happened to be checked by hand.

### R5 — Don't merge polylines to tidy geometry; merge only with a stated reason

Merging erases structure. The step-over density of a fault trace is a
physical quantity — Qiao et al. (2022) use it as a maturity proxy — so
concatenating polylines destroys the very thing being modelled.

Legitimate reasons, both used for Xianshuihe: the junction is far below the
step-over cut-off (a digitisation break, not a structure), or a fault is too
short to resolve. Record the reason in the compset README.

### R6 — Parallel strands cannot be merged into a chain

End-to-end joins collapse into one fault. Side-by-side strands do not — they
need multi-surface topology (`subei.gmsh.lite`), not chain topology
(`gulang.gmsh.lite`).

### R19 — Every compset parameter must be assigned to `par`

`case.setup` reads attributes off the `par` object, so a bare
`name = value` in `user_defined_params.py` sets a module-level local that
nothing ever reads.

**Why:** `fric_fini = 0.45` was missing the prefix in `gulang`, `subei`
and `xianshuihe`. It never had any effect. It went unnoticed only because
the value happens to equal the default in `defaultParameters.py` -- change
the number and nothing would happen, silently.

---

## Loading conventions

### R7 — `ftLoadMaxShear` is a strain rate that sets a stress, not a slip rate

`interstress.f90` computes `rs = rd·cos(2θ)·ant`, an asymptotic **stress**.
The long-term slip rate is emergent, so a published slip rate cannot be
substituted directly, and inverting `γ̇ = V/2W` against the mesh width gives
nonsense (~570 mm/yr for the SAF).

Calibrate: scale from the SAF anchor (`1.427e-14 s⁻¹` → ~38 mm/yr peak),
run, measure with `plot_saf_figure3.py`, rescale linearly.

### R8 — `ftLoadAngle = -999` is a placeholder, not a measurement

It makes the code assume loading is along x and derive θ from the fault
tangent. Fine where the rotated frame aligns with plate motion; elsewhere it
discards the real orientation of the strain field. A proper strain-rate map
supplies both the max shear magnitude and its direction.

### R9 — Domain extension and `term` are coupled

`ext = 10 + totalSimuTime*vp/1e3` km keeps boundary reflections from reaching
the fault before the rupture ends. Shrinking `ext` for a cheaper mesh breaks
that guarantee — at `ext = 100 km` a P-wave round trip is ~33 s against a
200 s run. Acceptable for meshing and geometry checks; state it in the file
and restore before trusting rupture dynamics.

---

## Reproducing published results

### R10 — Parity is against the reference implementation running, not its printed output

Gate a port on executing the original (MATLAB is on PATH) against the
published data. Printed tables are rounded and can be wrong.

**Why:** published Table 2 prints 49 events at Frazier Mountain where its own
script on its own published data computes 52. Matching the paper would have
meant breaking parity with the code that made it.

### R11 — Reproduce reference bugs deliberately, and say so

MATLAB's `load()` reads gfortran's dropped-E `0.8384675-101` as `0.8384675`
— truncation, not scientific notation. Reproducing published numbers requires
truncating too. Document any such choice at the point of use; never silently
"fix" a reference implementation while claiming parity.

### R12 — State the magnitude convention wherever magnitudes are compared

`plotRuptureDynamics` uses the model's own `mu = rou*vs^2 = 3.204e10 Pa`; the
paper's Figures 6 and 9 use `3500^2*3000 = 3.675e10 Pa`. Every catalog
magnitude and b-value therefore sits 0.04 Mw below the paper's scale.

### R13 — Constants derived from an external chain must record their derivation

`OBSERVED_EQDYNA_X_KM` was wrong by up to 147 km, misassigning a third of the
SSAF nodes to the wrong observed-rate band in every Figure 3 ever produced.
It carried no comment saying where it came from, so nothing flagged it.
Derive, verify against the source chain, and record how.

### R14 — Hand-picked cycle ids do not transfer between runs

Sequences are chaotic. Selecting events by cycle id from another run
mislabels unrelated events — Figure 9 stamped "1857-like" on an M5.7 that
broke 8 km of fault. Select by rupture footprint and label only on a match.

---

## Case I/O

### R15 — Read case outputs from the case dir **or** `aRawSimuData/`

`run.sh` moves `totalop.txt*`, `cyclelog.txt*`, `interval.txt*` (and
`binaryop` for C_mesh=2) into `aRawSimuData/` when a run finishes. A script
that only looks in the case dir breaks on every completed run.

Search both, and never sum both — double-counting `interval.txt` inflated the
simulated duration.

### R16 — `run.sh` must not delete `binaryop` before launching

It is the restart state. Deleting it made the documented
restart-from-binaryop procedure impossible for C_mesh=2 cases.

### R17 — Detach the run wrapper

`setsid` where available, so a process-group kill cannot orphan the binary
and skip the post-run steps. A `paper.saf.A` run lost its wrapper mid-flight:
the binary finished all 4000 cycles, but nothing was plotted, nothing was
moved to `aRawSimuData/`, and no `Job finished` marker was written.
