# Release Log v2.0.4

EQdyna.2Dcycle Release Documentation

This document details all changes and improvements in version 2.0.4.

## [2.0.4] - 2025-02-04

### Added
- Comprehensive workflow documentation in README.md
- MIT License (LICENSE file)
- Complete Quick Start guide with validated commands
- Enhanced installation script with macOS environment support
- Automated post-processing with plot generation
- Git ignore configuration for work directories
- Proper citations for academic references

### Enhanced
- Installation script (`install.eqdyna.2dcycle.sh`) with macOS environment setup
- Case creation workflow with environment variable support
- Python script compatibility and path handling
- Mesh generation library (`meshGenLib.py`) improvements
- Default parameters configuration
- Plot generation for rupture dynamics

### Fixed
- Environment variable handling in case creation scripts
- Python path resolution for case setup
- File permissions for executable scripts
- Mesh generation stability

### Documentation
- Complete README.md with step-by-step workflow
- Verified academic citations:
  - Duan & Oglesby (2006) - Core methodology
  - Liu et al. (2022) - San Andreas fault applications
- Project structure documentation
- Dependencies and requirements

### Validated
- Complete workflow from installation to post-processing
- Multi-cycle earthquake simulation capability
- Test case: Subei fault system (10 earthquake cycles)
- Magnitude range: M5.6-M7.3+ earthquakes
- Inter-event periods: 1-618 years
- Automatic visualization generation

### Technical Details
- **Fortran Compiler**: gfortran
- **Python Dependencies**: numpy, matplotlib, xarray
- **Optional Tools**: gmsh, meshio, nbconvert
- **Simulation Engine**: 2D Finite Element Method
- **Friction Laws**: Multiple supported friction models
- **Output Formats**: Text files, PNG plots

### Example Results
- Successfully simulates realistic earthquake cycles
- Stress evolution and rupture dynamics
- Fault slip and rupture time analysis
- Magnitude-frequency relationships

## Development Notes

### Repository Structure
- `src/` - Core Fortran 90 simulation engine
- `scripts/` - Python utilities for case management
- `compset/` - Pre-configured test cases
- `work/` - User simulation workspace (git-ignored)

### Supported Platforms
- macOS (tested with Homebrew)
- Ubuntu/Linux
- Manual compilation available for other Unix-like systems

---

**Note**: This changelog represents the current stable release. For development history, see git commit log.