# paper.saf.A

Paper-faithful reproduction of Liu et al. (2022) JGR Solid Earth Model A —
Southern San Andreas + San Jacinto fault system. C_mesh=2 (internal Fortran
structured-quad mesh, `meshgen1.f90` builds it from the per-fault
`x{1,2,3}_1.txt` geometry files shipped here — no `meshgen.py` step). Standard
build/run flow: repo root `README.md` § Run (`create.newcase` → `case.setup` →
`bash run.sh`; skip the `meshgen.py` step). `OMP_NUM_THREADS` defaults to 1;
default `icend = 4000` cycles — Ctrl-C / `kill` whenever you have enough for
the analysis at hand.

## Loading

`Rate_direction.txt`: per-fault-node loading, column 1 = loading rate (s^-1),
column 2 = angle between loading direction and fault tangent (deg).

**Provenance caveat.** This is the smooth per-node table recovered from a
paper-era **local** archive. It is not the file in either public archive, and
it is not confirmed to be the one that produced the published Model A results.
Neither the Zenodo nor the Pangaea copy of `Rate_direction.txt` can be that
file either — both are in an nrad/yr literal convention that the compiled
binary's `interstress.f90` would read as ~1e24 Pa. See
[`PROVENANCE.md`](PROVENANCE.md) for the full runtime trace. The authoritative
loading file for the published run appears not to have been archived anywhere.

## Files

```
user_defined_params.py    # case parameters (icend, friction, viscosity, ...)
x1_1.txt .. x3_1.txt      # per-fault structured-mesh geometry
Rate_direction.txt        # per-node loading rate + angle
RATE_DIRECTION.md         # that file's format, units and code path
PROVENANCE.md             # published-archive trace behind the caveat above
README.md
```

## Validation against the published results

Post-processing is reproduced by the Python ports in `scripts/`, each gated on
running the published MATLAB (Zenodo 10.5281/zenodo.5823021) headless against
the published Pangaea output (10.1594/PANGAEA.940262). Fetch both with `bash
scripts/fetch_published_reference.sh`; render the whole set for a case with
`python3 scripts/make_paper_figures.py <case_dir>`.

### 4000-cycle run here vs. published Model A

Rupture behaviour reproduces well; the timing does not.

| | this compset | published Model A |
|---|---|---|
| events / simulated time | 4000 / 10.76 kyr | 4000 / 15.17 kyr |
| Mmax | 8.23 | 8.27 |
| events M>7.0 | 88 (2.2%) | 92 (2.3%) |
| 1857-extent ruptures | 46 (1.1%) | 48 (1.2%) |
| SSAF-north-only ruptures | 314 (7.8%) | 382 (9.6%) |
| recurrence at FM | 225 yr, COV 0.24 | 283 yr, COV 0.41 |
| recurrence at WW | 235 yr, COV 0.27 | 308 yr, COV 0.42 |

Magnitude distribution and rupture-extent mix match; the interseismic clock
does not — the same 4000 events occur ~30% faster, and FM/WW are more
periodic. That points at the loading (consistent with the `Rate_direction.txt`
provenance caveat above), not at the dynamic rupture solver.

### Conventions that differ from the paper — read before comparing numbers

- **Rigidity in magnitude.** The paper's Figures 6 and 9 compute moment with `mu = 3500^2 * 3000 = 3.675e10 Pa` and a 22 km
  locking depth. `scripts/plotRuptureDynamics` uses `mu = 2670 * 3464^2 = 3.204e10 Pa`, the model's own density and shear-wave
  speed, at the same 22 km. Every magnitude and b-value it reports therefore sits `2/3*log10(3.675/3.204) = 0.04` Mw **below**
  the paper's published scale. `plot_saf_figure6.py` and `plot_saf_figure9.py` use the paper's constants, so their magnitudes
  are on the paper's scale and are NOT directly comparable to `catalog.csv`. The paper's constant is itself inconsistent with
  its own model parameters; making this selectable is planned.
- **Figure 4 windows.** The paper's Fig 4 panels start at 3 kyr and skip 0-3 kyr as burn-in (`tstart = [3, 6, 9, 12]`).
- **Figure 9 event selection.** The published MATLAB hard-codes cycle ids (856, 1041, 355, ...) hand-picked from the published
  sequence; those mean nothing in another run since sequences are chaotic, so `plot_saf_figure9.py` selects by rupture footprint
  instead and labels an event "1857-like" only when its extent actually matches. Use `--published` to reproduce the paper's
  exact panels on the archived data.

### Known error in the published Table 2

At Frazier Mountain the paper prints **49** counted events; the paper's own
script, run on the paper's own published data, computes **52**. The other five
FM statistics reproduce exactly (mean 283.43, std 115.64, COV 0.408, slip
8.697/3.268), and BF (107) and WW (48) match exactly. Node choice was swept
across the full fault-3 range and threshold sensitivity ruled out (all 52
events clear 0.5 m; the smallest is 0.74 m). `scripts/paleo_site_stats.py`
reports 52 and flags the delta rather than tuning to the printed value.

### Two bugs in the published post-processing, reproduced deliberately

- MATLAB's `load()` reads gfortran's dropped-E 3-digit exponents (`0.8384675-101`) as `0.8384675` — it truncates rather than
  reading scientific notation, verified against R2024b. The published statistics therefore treat those near-zero slips as ~0.84
  m. The ports truncate too, since that is what reproduces the published numbers; only paper-era output is affected, current
  runs emit no such tokens.
- `Figure5_Plot_Recurrene_Stats.m` never resets its moment accumulator per cycle and sums only faults 1-2, never fault 3 where
  all three paleo sites sit. Its magnitude output is unusable, but Table 2 does not consume it, so the port omits magnitude
  entirely.

## Companion compset

`compset/saf.gmsh.lite/` is the C_mesh=3 (gmsh) reproduction of the same
physics on a coarser unstructured mesh, useful for quicker exploratory runs.
