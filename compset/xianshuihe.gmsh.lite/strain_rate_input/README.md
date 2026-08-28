# Strain-rate input

`GSRM_v2.1_xianshuihe_region.txt` — regional subset of the Global Strain
Rate Model v2.1 (Kreemer, Blewitt & Klein, 2014, G-cubed 15, 3849-3889,
doi:10.1002/2014GC005407), cut to 99-104 E / 27.5-33 N on the model's
0.1 deg grid. 2750 cells.

**Not vendored** — the data is public and citable, so fetch it:

```bash
bash fetch_strain_rate.sh      # downloads, decompresses, cuts to the region
python3 strain_rate_loading.py # then computes the loading columns
```

The global file is 88 MB compressed / 529 MB expanded, from
<https://geodesy.unr.edu/GSRM/>; the script keeps only the ~240 KB
regional cut.

Columns, as documented in the GSRM header. **Units are 1e-9/yr**:

| col | field | note |
|---|---|---|
| 1 | lat | deg |
| 2 | lon | deg |
| 3 | exx | 1e-9/yr |
| 4 | eyy | 1e-9/yr |
| 5 | exy | 1e-9/yr |
| 6 | vorticity | |
| 7 | RL-NLC | right-lateral, no length change |
| 8 | LL-NLC | left-lateral, no length change |
| 9 | e1 | principal, 1e-9/yr |
| 10 | e2 | principal, 1e-9/yr |
| 11 | azi_e1 | azimuth of e1, deg |

Licence: CC-BY-NC-SA 3.0 (unported), GEM Foundation. Cite Kreemer et al.
(2014) on any use.

## Status

`strain_rate_loading.py` currently samples the **123 KML control points**,
not mesh fault nodes (the 5-fault mesh has 1026). The output CSV is a
profile of the field, **not yet `nsmpGeoPhys.txt` input**. Two things are
needed before it is:

1. Resample onto the mesh fault nodes (`nsmp.txt` + `vert.txt`). The frames
   already match — both use `principal_frame()` from `plot_fault_trace.py`.
2. Interpolate rather than take the nearest GSRM cell. The grid is 0.1 deg
   (~9-11 km); at 400 m node spacing ~25 consecutive nodes would otherwise
   share one cell value.
