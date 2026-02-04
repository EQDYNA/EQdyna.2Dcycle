# EQdyna.2Dcycle

EQdyna.2Dcycle is a 2-D finite element method (FEM) code to simulate physics-based multicycle earthquake dynamics on geometrically complex fault systems.

## Quick Start

### 1. Install and Build

For macOS:
```bash
./install.eqdyna.2dcycle.sh -e macos
```

For Ubuntu/Linux:
```bash
./install.eqdyna.2dcycle.sh -m ubuntu
```

Manual build:
```bash
cd src/
make
cd ..
mkdir -p bin
mv src/eqdyna.2dcycle bin/
```

### 2. Set Environment Variables

```bash
export EQDYNA2DCYCLEROOT=$(pwd)
export PATH="$EQDYNA2DCYCLEROOT/bin:$EQDYNA2DCYCLEROOT/scripts:$PATH"
```

### 3. Create and Run a Test Case

```bash
# Create new case
mkdir -p work
python3 scripts/create.newcase work/test_subei test.subei

# Configure case
cd work/test_subei
python3 case.setup

# Run simulation
chmod +x run.sh
./run.sh
```

### 4. View Results

Results are automatically generated in the `aPlots/` directory, including rupture dynamics visualizations.

## Example Output

A successful run simulates multiple earthquake cycles:
- Major earthquake (M7.3+) 
- Multiple smaller earthquakes (M5.6-5.8)
- Inter-event times ranging from 1-600+ years
- Automatic post-processing with plots

## Project Structure

- `src/` - Fortran 90 source code
- `scripts/` - Python utilities for case management
- `compset/` - Pre-configured test cases
- `work/` - User simulation cases (git-ignored)

## Dependencies

- **Fortran**: gfortran compiler
- **Python 3**: numpy, matplotlib, xarray
- **Optional**: gmsh, meshio, nbconvert

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citations

If you use EQdyna.2Dcycle in your research, please cite:

**Core methodology:**
- Duan, B., & Oglesby, D. D. (2006). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike‐slip faults. *Journal of Geophysical Research*, 111(B5), B05309. https://doi.org/10.1029/2005JB004138

**San Andreas fault applications:**
- Liu, D., Duan, B., Scharer, K., & Yule, D. (2022). Observation‐constrained multicycle dynamic models of the southern San Andreas and the northern San Jacinto faults: Addressing complexity in paleoearthquake extent and recurrence with realistic 2D fault geometry. *Journal of Geophysical Research: Solid Earth*, 127(2), e2021JB023420. https://doi.org/10.1029/2021JB023420

## Contact

- Website: https://seismotamu.wixsite.com/emlam
- Contacts: dliu@ig.utexas.edu, bduan@tamu.edu
- Product of Earthquake Modeling Lab @ Texas A&M University
