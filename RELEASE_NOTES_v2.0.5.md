# EQdyna.2Dcycle v2.0.5 Release Notes

**Release Date**: February 5, 2025

## 🎉 Major Improvements

### 🔄 Complete Testing System Overhaul
- **Automated Test Pipeline**: `python3 -m test.testAll` runs complete end-to-end testing
- **Integrated Verification**: Automatic comparison against reference results
- **Modular Architecture**: Separate verification script (`test/verify.test.py`)
- **SUCCESS/FAIL Reporting**: Clear validation status for all output files

### 📁 Organized Project Structure
- **Clean Folder Organization**: 
  - `user_fault_geometry_input/` for input fault geometries
  - `fem_mesh_output/` for generated mesh files
  - `aRawSimuData/` for simulation results
- **Smart File Management**: Automated copying of mesh files during simulation
- **Git Integration**: Added `bin/` to `.gitignore` to prevent executable commits

### 🏗️ Enhanced Case Management
- **Improved create.newcase**: Handles organized directory structures
- **Smart run.sh Generation**: Automatically copies required mesh files
- **Streamlined Workflow**: From fault geometry to results in organized folders

### 🧪 Robust Meshing Workflow
- **GMSH Integration**: High-quality quadrilateral mesh generation
- **Fault Interface Support**: Proper split-node implementation for fault mechanics
- **Organized Output**: Mesh files cleanly separated in `fem_mesh_output/`
- **Automatic Validation**: Required file checking before simulation

## 🛠️ Technical Enhancements

### Testing Framework
```bash
# Complete automated testing
python3 -m test.testAll

# Manual verification
python3 test/verify.test.py
```

**Test Coverage:**
- ✅ Fresh compilation from source
- ✅ Case creation with organized structure
- ✅ Mesh generation and validation
- ✅ Full 10-cycle earthquake simulation
- ✅ Post-processing and visualization  
- ✅ Verification against reference data

### Improved Workflow
```bash
# New streamlined process
create.newcase work/my_case test.subei
cd work/my_case
python3 meshgen.py          # Creates fem_mesh_output/
python3 case.setup          # Generates smart run.sh
bash run.sh                 # Runs simulation with automatic file management
```

### Code Organization
- **Renamed Components**: `test.subei.dev` → `test.subei` for clarity
- **Flexible Imports**: Support for both module and standalone execution
- **Error Handling**: Comprehensive validation and user feedback
- **Path Management**: Automatic project root detection

## 📊 Validation Results

### Test Verification Status
All reference comparisons: **SUCCESS**
- ✅ `cyclelog.txt1` - Earthquake cycle information
- ✅ `interval.txt1` - Inter-event intervals  
- ✅ `totalop.txt1-10` - Complete stress/slip time series

### Performance Metrics
- **Build Time**: ~15 seconds (fresh compilation)
- **Mesh Generation**: ~2 seconds (13,538 nodes, 14,236 elements)
- **Simulation Runtime**: ~45 seconds (10 earthquake cycles)
- **Total Testing Time**: ~67 seconds (complete pipeline)

## 🔧 System Requirements

### Dependencies
- **Fortran**: gfortran compiler
- **Python 3**: numpy, matplotlib, xarray
- **Mesh Generation**: gmsh (for advanced geometries)
- **Testing**: Python module support

### Platform Support
- ✅ macOS (tested)
- ✅ Linux/Ubuntu  
- ✅ Unix-like systems

## 📝 File Structure Changes

### New Organization
```
EQdyna.2Dcycle/
├── test/                           # Testing framework
│   ├── testAll.py                 # Main test runner
│   ├── verify.test.py             # Result verification
│   ├── testNameList.py            # Test configuration
│   └── reference.results/         # Reference data
├── compset/test.subei/            # Renamed from test.subei.dev
└── work/                          # Simulation workspace
    └── case_name/
        ├── user_fault_geometry_input/  # Input geometries
        ├── fem_mesh_output/           # Generated mesh
        └── aRawSimuData/              # Simulation results
```

### Updated Scripts
- **case.setup**: Generates smart `run.sh` with automatic file copying
- **create.newcase**: Handles organized directory structures
- **meshgen.py**: Outputs to organized `fem_mesh_output/` folder

## 🚀 Migration Guide

### For Existing Users
1. **Update Testing**: Replace `python3 check.test.py` with `python3 -m test.testAll`
2. **Organized Structure**: New cases automatically use organized folders
3. **Mesh Files**: Now cleanly separated in `fem_mesh_output/`
4. **Git Workflow**: `bin/` folder no longer tracked

### For New Users
Follow the streamlined workflow in README.md - everything is automated!

## 📖 Documentation Updates

### Enhanced CLAUDE.md
- Updated testing instructions
- Corrected workflow commands
- Added organized folder descriptions

### README.md Updates
- New automated testing section
- Organized workflow documentation
- Updated command examples

## 🎯 Key Benefits

1. **One-Command Testing**: Complete validation with `python3 -m test.testAll`
2. **Clean Organization**: Logical separation of inputs, mesh, and outputs
3. **Automated Workflows**: Smart file management reduces manual steps
4. **Robust Validation**: Comprehensive verification against reference data
5. **Developer Friendly**: Modular testing architecture for easy extension

## 📧 Support

- **Website**: https://seismotamu.wixsite.com/emlam
- **Contact**: dliu@ig.utexas.edu, bduan@tamu.edu
- **Institution**: Earthquake Modeling Lab @ Texas A&M University

---

**License**: MIT License  
**Platform Support**: macOS, Linux, Unix-like systems  
**Python Compatibility**: 3.7+