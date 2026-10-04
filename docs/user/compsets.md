# Compsets

A compset is a ready-made case: fault geometry, loading, parameters and, for
C_mesh = 3, a mesh generator. Each has its own README in
`compset/<name>/` with the details.

| compset | mesh | fault system |
|---|---|---|
| `paper.saf.A` | C_mesh = 2 | southern San Andreas + northern San Jacinto, Liu et al. (2022) Model A |
| `saf.gmsh.lite` | C_mesh = 3 | the same system on a coarser unstructured mesh, for quicker runs |
| `subei.gmsh.lite` | C_mesh = 3 | Subei fault system (atf, dxs, sbt) |
| `gulang.gmsh.lite` | C_mesh = 3 | Gulang fault system, 5 faults |
| `xianshuihe.gmsh.lite` | C_mesh = 3 | Xianshuihe fault, eastern Tibet, 7 faults, loading from GSRM v2.1 |

## The two meshing modes

| | C_mesh = 2 | C_mesh = 3 |
|---|---|---|
| mesh | structured quads, built by the solver | unstructured quads, built by gmsh |
| geometry input | `x*_1.txt` per fault | `user_fault_geometry_input/*.gmt.txt` |
| loading input | `Rate_direction.txt` | `nsmpGeoPhys.txt`, written by `meshgen.py` |
| extra step | none | `python3 meshgen.py` |

Both are first-class: at a fixed thread count, `paper.saf.A` and
`saf.gmsh.lite` each reproduce their stored reference results exactly.

## Notes per compset

**`paper.saf.A`.** Reproduces the published Model A rupture behaviour; its
recurrence clock runs about 30% faster than the published sequence. The
loading file's provenance is documented in
`compset/paper.saf.A/PROVENANCE.md`.

**`xianshuihe.gmsh.lite`.** Loading is sampled from the GSRM v2.1 strain-rate
model at the mesh fault nodes. The data are not shipped; fetch them with
`bash compset/xianshuihe.gmsh.lite/fetch_strain_rate.sh`. Use a target
asymptotic shear stress of 90–100 MPa (see [Model inputs](inputs.md)).
