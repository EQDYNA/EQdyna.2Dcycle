# San Andreas Fault System Meshing (Development Version)

This directory contains the development version of the meshing system for the Southern San Andreas Fault (SAF) system used in EQdyna.2Dcycle earthquake dynamics simulations.

## Overview

The meshing system generates 2D finite element meshes for the SAF fault system using GMSH, then processes them for use with the EQdyna.2Dcycle simulation engine. It handles fault geometry, creates split nodes for fault interfaces, and assigns physics parameters.

## Fault System Description

The **SAF system** consists of three disconnected faults:
- **ssaf**: Southern San Andreas Fault (80 control points, main fault trace)
- **sjfn**: San Jacinto Fault North Strand (15 control points)
- **sjfs**: San Jacinto Fault South Strand (10 control points)

All three are right-lateral strike-slip faults with vertical dip (90 degrees).

For mesh topology purposes, SSAF is split into two segments (ssaf1: 67 pts, ssaf2: 13 pts) at the branch point nearest to SJFN. In the simulation, ssaf1 and ssaf2 are recombined as one fault 'ssaf'.

### Changes from test.subei

| Component | test.subei | test.saf.dev |
|---|---|---|
| Fault system | Subei (atf1, atf2, dxs, sbt) | SAF (ssaf, sjfn, sjfs) |
| Fault types | Mixed (strike-slip + thrust) | All right-lateral strike-slip |
| ftNamesForGmsh | ['atf1', 'atf2', 'dxs', 'sbt'] | ['ssaf1', 'ssaf2', 'sjfn', 'sjfs'] |
| ftNames | ['atf', 'dxs', 'sbt'] | ['ssaf', 'sjfn', 'sjfs'] |
| system | "subei" or "none" | "saf" |
| dx (mesh size) | 0.5 km | 0.3 km |
| totalSimuTime | 15 s | 30 s |
| Supplementary point | supPt1 used | Not needed |
| Fault connectivity | Faults connected at junctions | Faults disconnected |

### Mesh Topology (6 surfaces)

1. **Top block**: Above SSAF (ssaf1 + ssaf2 + boundary)
2. **Left block**: Triangle left of ssaf1 start
3. **SSAF-SJFN wedge**: Between ssaf2 and sjfn
4. **Bottom-left block**: Below ssaf1, left of sjfn start
5. **SJF-bottom block**: Below sjfn and sjfs to bottom boundary
6. **Right block**: Right of ssaf2/sjfn/sjfs to right boundary

## Files and Structure

### Input Files
```
user_fault_geometry_input/
├── ssaf1.gmt.txt         # SSAF segment 1 (67 pts)
├── ssaf2.gmt.txt         # SSAF segment 2 (13 pts)
├── sjfn.gmt.txt          # San Jacinto Fault North Strand (15 pts)
└── sjfs.gmt.txt          # San Jacinto Fault South Strand (10 pts)

user_defined_params.py    # EQdyna simulation parameters
userDefinedFaultSysGeoPhys.py  # Fault physics definitions (system="saf")
```

### Meshing Tools
```
meshgen.py                # Python script for mesh generation
meshGenLib.py             # Meshing utility library
meshgen.ipynb             # Jupyter notebook (needs updating for SAF)
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
ftType = -1           # Right-lateral strike-slip (all 3 faults)
ftDip = 90            # Vertical dip
ftLoadMaxShear = 1.427e-14
ftVis = 6e21          # Pseudo viscosity (Pa-s)
```

## Dependencies

```bash
pip install numpy matplotlib meshio gmsh
```

## Fault Geometry Source

Fault trace coordinates were converted from the previous version's `x*.txt` files:
- `x3_1.txt` (80 pts) -> ssaf (Southern San Andreas Fault)
- `x1_1.txt` (15 pts) -> sjfn (San Jacinto Fault North Strand)
- `x2_1.txt` (10 pts) -> sjfs (San Jacinto Fault South Strand)

---

**Last Updated**: February 2025
