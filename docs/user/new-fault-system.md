# Setting up a new fault system

The steps that take a digitised fault trace to a running model. Start from
the compset nearest to yours (`xianshuihe.gmsh.lite` is the template for
loading from a strain-rate field) and copy its scripts; each has its region
and file names at the top.

## 1. Literature first

Write `compset/<name>/LITERATURE.md` before touching the loading. It fixes
what the model will be judged against, so the targets cannot drift towards
whatever the model produces. Record, each with its source:

- slip sense and slip rates per fault, geological and geodetic, and which
  parallel faults share the motion but are left out of the model;
- the historical and paleoseismic record: magnitudes, dates, intervals,
  rupture extents;
- seismogenic depth and coupling (locked segments, creeping barriers);
- earlier physics-based models of the same system;
- a closing table of targets: largest magnitude, slip per large event,
  recurrence, how often ruptures cross the step-overs.

## 2. Geometry

`export_fault_geometry.py` turns a lon/lat trace into `ft*.gmt.txt` in km:

- **Frame:** origin at the trace centroid, x along the principal axis of the
  points (SVD). meshgen fits y(x) per fault, so every fault must be monotone in
  x after rotation; the exporter fails if one is not.
- **Resample** each polyline to even arc length. Raw digitisations jitter by
  metres and can double back; the resampled points are what meshgen receives.
- **Merge** only sub-kilometre collinear breaks; keep real step-overs as
  separate faults.
- **Number** the faults in one direction along the system and say which in the
  README.

Set `ntotft` and `ftcn` (control points per fault) in `user_defined_params.py`
and add the system's branch to `userDefinedFaultSysGeoPhys.py` and `meshgen.py`.

## 3. Mesh

```bash
create.newcase --work_dir work/<case> --compset <name>
cd work/<case> && python3 case.setup && python3 meshgen.py
python3 checkMeshQuality.py . && python3 plotMeshFaults.py .
```

Accept the mesh when it has no triangles, no orphaned split nodes and clean
fault tips, with interior angles roughly 40–145° and aspect ratio below ~4.
Put the counts in the compset README.

## 4. Slip sense

Set `par.slipSense` in `user_defined_params.py`: +1 for right-lateral, −1 for
left-lateral. **`ftType` does not set it**; the solver never reads `ftType`.
One value covers the whole system.

## 5. Loading from a strain-rate field

1. **Fetch** the tensor field for a box around the faults:
   `bash fetch_strain_rate.sh` (GSRM v2.1; set the lat/lon box at the top).
   The data are third-party, so fetch them, never commit them; add the file to
   `.gitignore`.
2. **Sample at the mesh nodes:**
   `python3 strain_rate_loading.py --case work/<case>` writes
   `<name>_strain_loading.csv` (γ and φ per node) and the Figure-2-style
   panels. Interpolate the tensor components, never the principal angle.
3. **Audit** what it prints:
   - *angle vs resolved tensor* error near round-off: φ is defined so that
     γ cos 2φ is the resolved shear in the `slipSense` direction and
     γ sin 2φ the deviatoric normal strain rate;
   - *resolved shear in the slip direction* at nearly every node: if not, the
     slip sense or the field is wrong;
   - *largest T* that keeps every node compressive.

   Also check γ against the slip rates in `LITERATURE.md`: a locked fault
   slipping at V with locking depth D has peak surface shear strain rate about
   V / (2πD), so V ≈ 2πDγ. The field is smoothed over its grid cells, so expect
   this to come out at or a little below the geological rate.
4. **Choose the target stress T** (see [Model inputs](inputs.md#choosing-the-loading-stress)).
   The asymptotic stresses are shear T cos 2φ and normal −N − s·T sin 2φ. A
   node can nucleate only where T cos 2φ exceeds f_s times its normal stress;
   strongly clamped strands never nucleate and slip only when a rupture
   reaches them. Pick two or three values up to the printed limit.
5. **Apply:** `python3 apply_strain_loading.py --case work/<case> --target-stress <T in Pa>`
   writes γ, φ and η = T/γ into `fem_mesh_output/nsmpGeoPhys.txt`; `run.sh`
   copies it into the case.

## 6. Sweep, then run

Run one case per T value, with the same thread count throughout (results are
bit-reproducible only at a fixed `OMP_NUM_THREADS`). Two hundred cycles are
enough to see whether the system settles. Compare each against the
`LITERATURE.md` targets with the catalogue and figures
([Post-processing](post-processing.md)), choose T, then run the production
sequence.

## 7. Stage

Commit the compset's scripts, `LITERATURE.md`, the README (mesh counts,
loading audit, chosen T and why) and the three figures worth keeping: the
trace, the strain map and the on-fault loading. Leave out the fetched data and
the per-mesh CSV, which is stale as soon as the mesh changes. Add a row to
[Compsets](compsets.md).
