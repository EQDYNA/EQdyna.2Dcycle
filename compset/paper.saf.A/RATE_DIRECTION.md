# `Rate_direction.txt` — definition, units, and how it drives interseismic stress

Applies to `C_mesh = 2` cases, i.e. `compset/paper.saf.A`. The C_mesh=3
equivalent is the 9-column `nsmpGeoPhys.txt`; see `CLAUDE.md`. For where the
paper's own file came from and why its provenance is uncertain, see
[`PROVENANCE.md`](PROVENANCE.md).

This file is the per-node geodetic forcing for C_mesh=2 cases. It carries
the **maximum shear strain rate** γ(x) and the **angle of compression** φ(x)
at every fault node, which together define the asymptotic interseismic
loading via paper Liu et al. (2022) eqs (1)–(3).

## 1. Physical meaning

Per the Smith‑Konter & Sandwell (2009) regional strain‑rate field,
each on‑fault node sees:

- **γ** = maximum shear strain rate at that location (rad/s)
- **φ** = angle of compression = (local fault strike) − (max shear strain‑rate direction), measured in degrees

These are the inputs to paper eq (3):

```
γ_τ = γ cos(2φ)       (shear-strain-rate component on the fault)
γ_n = γ sin(2φ)       (normal-strain-rate component on the fault)
```

and then to eqs (1)–(2), which evolve fault stress over an interseismic period:

```
σ_τ(t) = (σ_τ⁰ − η γ_τ)·exp(−μ t / η) + η γ_τ                    (1)
σ_n(t) = (σ_n⁰ − σ^a − η γ_n)·exp(−μ t / η) + η γ_n + σ^a        (2)
```

where μ = shear modulus, η = viscosity, σ^a = ambient (lithostatic) normal
stress, and σ_τ⁰, σ_n⁰ are the on‑fault stresses at the start of the
interseismic period.

## 2. Units (as consumed by current `src/interstress.f90`)

| Column | Symbol | Unit | Notes |
|---|---|---|---|
| 1 | γ | **rad/s** | typical magnitudes 1e‑15 to 2e‑14 |
| 2 | φ | **degrees** | per‑node "angle of compression"; ranges roughly −22° to +30° on the SAF |

Convert nrad/yr to rad/s by multiplying by `1e‑9 / (365.25·86400) ≈ 3.169e‑17`.

> **Heads up:** an older code path (`strbld.f90` / `strini.f90`, never
> linked into the runtime binary today) used `ant = ant0·450/rd(1)` with
> `rd(1)` in **nrad/yr literal** (e.g. `160`, `470`). The current
> `interstress.f90` uses **rad/s** and the simpler `ant = ant0·str/rd(1)`.
> If you find a `Rate_direction.txt` with col1 values like `160` or `470`,
> it's nrad/yr — convert before feeding the current binary.

## 3. File layout

```
maxftnode = 1769     ! max fault-node count across faults
ntotft    = 3        ! number of faults
total rows = maxftnode * ntotft = 5307
```

Each fault gets a slot of `maxftnode` rows. Only the first `nfnode(ift)`
rows of each slot are *active* — the rest are padding (typically `0 0`).

```
rows      1 ..  295   active fault 1 (NSJF Claremont, nfnode(1)=295)
rows    296 .. 1769   padding for fault 1 slot              (0 0)
rows   1770 .. 1947   active fault 2 (NSJF Clark, nfnode(2)=178)
rows   1948 .. 3538   padding for fault 2 slot              (0 0)
rows   3539 .. 5307   active fault 3 (SAF, nfnode(3)=1769)
```

The binary reads it in `eqdyna2d.f90:65–66`:

```fortran
allocate(rd(2, maxftnode*ntotft))
read(3,*) (rd(1,i), rd(2,i), i=1, maxftnode*ntotft)
```

…and the on‑fault loop in `interstress.f90` accesses `rd(:, ida(i))` with

```fortran
ida(ntag) = (ift - 1)*maxftnode + i_local
```

so only active rows are touched at runtime.

## 4. Code mapping (`src/interstress.f90`, C_mesh=2 branch)

| Paper symbol | Code variable / expression |
|---|---|
| μ | `amu = vs**2 * rou` (≈ 3.20e10 Pa) |
| η | `ant`, computed per node |
| σ_τ⁰ | `ss0(i)` — for ic==1 set to `−ambientnorm·fric_fini` |
| σ_n⁰ | `ns0(i)` — for ic==1 set to `ambientnorm` |
| σ^a | `ambientnorm` (negative; e.g. −100 MPa) |
| γ(x) | `rd(1, ida(i))` — *per‑node* from Rate_direction.txt col 1 |
| φ(x) | `rd(2, ida(i))` — *per‑node* from col 2 |
| γ_τ | `rd(1)·cos(2θ)` (with `θ = rd(2)·π/180`, clamped to ≤ 45°) |
| γ_n | `rd(1)·sin(2θ)` |
| η γ_τ | `rs` |
| η γ_n | `−rn` (note sign convention) |

Concretely:

```fortran
theta = rd(2,ida(i))/180.0d0*pi
if (theta >= 45.0d0/180.0d0*pi) theta = 45.0d0/180.0d0*pi
if (rd(1,ida(i)) > 0.0d0) then
    ant = ant0 * str / rd(1,ida(i))      ! η inversely ∝ γ
else
    ant = ant0                           ! fallback at endpoints with γ=0
endif
rs =  rd(1,ida(i)) * cos(2.0d0*theta) * ant
rn = -rd(1,ida(i)) * sin(2.0d0*theta) * ant
shs(i) = (ss0(i) - rs)             * exp(-tinter*amu/ant) + rs
ns(i)  = (ns0(i) - ambientnorm - rn) * exp(-tinter*amu/ant) + rn + ambientnorm
```

After the cancellation `rd·(ant0·str/rd) = str·ant0`, the asymptotes are:

```
σ_τ(∞) = η γ_τ = str · cos(2φ) · ant0       ← independent of γ at fixed φ
σ_n(∞) = σ^a + η γ_n = ambientnorm − str · sin(2φ) · ant0
```

That's why **spatial variation in σ_τ along strike comes from per‑node
variation in φ** (cos(2φ)). The relaxation time `τ = ant/μ` *does*
depend on γ, so high‑γ nodes approach the asymptote faster than low‑γ
nodes — but the asymptote itself is a function of φ alone.

## 5. Why per‑node smoothness matters

If φ is uniform across many adjacent nodes (as in the simplified
"stepped" Pangaea archived file, which has only φ ∈ {0°, 5.2°, 13°}),
then σ_τ(∞) is identical across those nodes. They all approach the same
threshold at similar times → when one nucleates, all neighbors are at
threshold → **cascade rupture** (M ~7.7 in cycle 1 instead of paper's
M ~5.8).

If φ varies smoothly per node (as in the paper's actual input — see
`source_local_recovered/input/Rate_direction.txt`, which has 295 unique
φ values on F1 and 1769 on F3), σ_τ(∞) varies smoothly along strike →
neighbors of the nucleation point are *not* at threshold → **localized
rupture** (M ~5.7).

We verified this end‑to‑end: with the smooth per‑node file, omp4
cycle 1 produces M = 5.74 with nucleation at node 475, slip 0.26 m
spanning ~5 km — matching Pangaea Model A cycle 1 (M5.76, 0.26 m) to
within 0.03 magnitude units.

## 6. Provenance of the two `Rate_direction.txt` files in this repo

| Path | Format | What it is |
|---|---|---|
| `docs/published/source_local_recovered/input/Rate_direction.txt` | **rad/s**, smooth per‑node, varying φ | The paper's *actual* runtime input. Recovered from local disk by D. Liu after publication. |
| `docs/published/pangaea_results/mesh/Rate_direction.txt` | **nrad/yr** literal, stepped chunks | Archived simplification on Pangaea (and identically in the Zenodo `dunyuliu/.../input/`). NOT what produced the published Pangaea outputs; appears to be a later quantized template. |
| `compset/paper.saf.A/Rate_direction.txt` | rad/s smooth (= source_local_recovered) | Working compset for paper‑A reproduction. |

To regenerate Fig 2 of the paper from a `Rate_direction.txt` plus the mesh,
see `work/paper_figure2_recovered.png` and the script that produced it
(`/usr/bin/python3` snippet, also in this doc's git history if needed).

## 7. Related runtime parameters from FE_Global.txt

```
ant0 = 8.4e21 Pa·s        ! reference viscosity (line 22). Combined with str/γ_max
                          !   gives the per-model η_min in paper Table 1.
str  = 1.427e-14 rad/s    ! "reference max shear strain rate" (line 23).
                          !   Equivalent to 450 nrad/yr; appears in η = ant0·str/γ.
ambientnorm = -100e6 Pa   ! σ^a (line 25)
fric_fini   = 0.45         ! initial fault shear/|normal| ratio (Model A; line 12)
                          !   sets σ_τ⁰ = −ambientnorm·fric_fini = 45 MPa
fric_fs     = 0.5          ! static friction (Model A; line 9). Determines threshold.
```

## 8. Common gotchas

- **Don't feed nrad/yr literal** (e.g. `160 13`) to the current binary. It
  treats col1 as rad/s, so 160 rad/s → asymptotic `rs = 160·8.4e21 ≈ 1.3e24` Pa
  → instant nucleation. Convert first (× 3.169e‑17) or use the smooth file.
- **The `0 0` padding rows are required** for the file to be `maxftnode·ntotft`
  long. The binary reads all of them but only consumes rows in active slots.
- **`nfnode(ift)` is set in `meshgen1.f90`** (C_mesh=2). If you regenerate the
  mesh and `nfnode` changes, you must regenerate `Rate_direction.txt` to
  match — the active‑row counts must agree.
