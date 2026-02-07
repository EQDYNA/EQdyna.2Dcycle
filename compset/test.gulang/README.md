# Gulang-Zhongwei Fault System Meshing

This directory contains the meshing system for the Gulang-Zhongwei fault system used in EQdyna.2Dcycle earthquake dynamics simulations.

## Overview

The meshing system generates 2D finite element meshes for a 5-fault system using GMSH, then processes them for use with the EQdyna.2Dcycle simulation engine. It handles fault geometry, creates split nodes for fault interfaces, and assigns physics parameters.

## Fault System Description

The **Gulang-Zhongwei system** consists of five nearly end-to-end left-lateral strike-slip faults forming a chain from west to east:
- **ft1**: 79 control points (westernmost)
- **ft2**: 53 control points (truncated from 63, see notes below)
- **ft3**: 83 control points
- **ft4**: 167 control points (longest)
- **ft5**: 34 control points (easternmost)

All five are left-lateral strike-slip faults with vertical dip (90 degrees).

### Changes from test.subei

| Component | test.subei | test.gulang |
|---|---|---|
| Fault system | Subei (atf1, atf2, dxs, sbt) | Gulang-Zhongwei (ft1-ft5) |
| Number of faults | 3 (simulation) / 4 (gmsh) | 5 |
| Fault types | Mixed (strike-slip + thrust) | All left-lateral strike-slip |
| ftNamesForGmsh | ['atf1', 'atf2', 'dxs', 'sbt'] | ['ft1', 'ft2', 'ft3', 'ft4', 'ft5'] |
| ftNames | ['atf', 'dxs', 'sbt'] | ['ft1', 'ft2', 'ft3', 'ft4', 'ft5'] |
| system | "subei" or "none" | "gulang" |
| dx (mesh size) | 0.5 km | 0.3 km |
| totalSimuTime | 15 s | 30 s |
| Supplementary point | supPt1 used | Not needed |
| Fault connectivity | Faults connected at junctions | Faults nearly end-to-end |

### Mesh Topology (4 surfaces)

The 5 faults form a nearly-collinear chain connected by short auxiliary lines. This gives a simple 4-surface topology:

1. **Top block**: Above the fault chain (chain + boundary top)
2. **Left triangle**: Left of ft1 start (ft1 start to left-top and left-bottom)
3. **Bottom block**: Below the fault chain (reversed chain + boundary bottom)
4. **Right triangle**: Right of ft5 end (ft5 end to right-top and right-bottom)

### Changes to meshGenLib.py

- **Array size fix**: Changed hardcoded `maxNumOfFtNodes*3` to `maxNumOfFtNodes*nFt` in `writeFilesForEQdyna()` to support arbitrary number of faults (lines 323-325).

### ft2 Truncation

The original ft2 (63 points) extended east past ft3's starting x-coordinate (ft2 end: x=-4.56, ft3 start: x=-5.51). This caused the auxiliary line F2_F3 (connecting ft2's end to ft3's start) to backtrack westward and cross ft2's own geometry, producing gmsh intersection errors.

**Fix**: ft2 was truncated to 53 points, cutting at the last point west of ft3's starting x-coordinate (point 53: x=-5.578). The auxiliary line from the truncated ft2 end to ft3's start is now a short 0.07 km segment that does not cross any fault.

## Coordinate Conversion

Fault traces were originally in lon/lat coordinates (shapefiles 1.shp through 5.shp). The script `convert_shp_to_gmt.py` converts them to local Cartesian km coordinates using:
- Reference point: lon=104.122925, lat=37.415949
- 1 deg lat = 111.32 km
- 1 deg lon = 111.32 * cos(lat) km

## Mesh Statistics

- **21,487 nodes**, **20,704 quad elements**
- Fault nodes: ft1 (207), ft2 (75), ft3 (111), ft4 (192), ft5 (133)
- Domain: ~485 km x ~388 km

## Files and Structure

### Input Files
```
user_fault_geometry_input/
├── ft1.gmt.txt           # Fault 1 (79 pts, westernmost)
├── ft2.gmt.txt           # Fault 2 (53 pts, truncated)
├── ft3.gmt.txt           # Fault 3 (83 pts)
├── ft4.gmt.txt           # Fault 4 (167 pts, longest)
└── ft5.gmt.txt           # Fault 5 (34 pts, easternmost)

convert_shp_to_gmt.py        # Shapefile to km coordinate converter
user_defined_params.py        # EQdyna simulation parameters
userDefinedFaultSysGeoPhys.py # Fault physics definitions (system="gulang")
```

### Meshing Tools
```
meshgen.py                # Python script for mesh generation
meshGenLib.py             # Meshing utility library (modified for 5 faults)
```

### Output Files
```
fem_mesh_output/
├── vert.txt              # Node coordinates - USED BY FORTRAN
├── fac.txt               # Element connectivity - USED BY FORTRAN
├── nsmp.txt              # Fault node mapping - USED BY FORTRAN
├── nsmpGeoPhys.txt       # Fault physics parameters - USED BY FORTRAN
├── meshGeneralInfo.txt   # Mesh statistics - USED BY FORTRAN
├── nsmpTanLen.txt        # Fault tangent and length data
├── eqdynaMesh.msh        # GMSH mesh file
├── meshWOSplitNode.png   # Mesh visualization
└── meshWithFaultNodes.png # Mesh with fault nodes highlighted
```

## Usage

```bash
# Convert shapefiles to .gmt.txt (only needed once)
python3 convert_shp_to_gmt.py

# Ensure fem_mesh_output/ directory exists
mkdir -p fem_mesh_output

# Run meshing pipeline
python3 meshgen.py
```

## Key Parameters

### Mesh Resolution
```python
dx = 0.3              # Fault zone mesh size (km)
dxAtBoundary = 20     # Boundary mesh size (km)
totalSimuTime = 30    # Simulation time for boundary sizing (s)
```

### Fault Physics (userDefinedFaultSysGeoPhys.py)
```python
ftType = 1            # Left-lateral strike-slip (all 5 faults)
ftDip = 90            # Vertical dip
ftLoadMaxShear = 1.427e-14
ftVis = 6e21          # Pseudo viscosity (Pa-s)
```

## Dependencies

```bash
pip install numpy matplotlib meshio gmsh geopandas shapely
```

---

**Last Updated**: February 2026
