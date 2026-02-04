# EQdyna.2Dcycle v2.0.4 Release Notes

**Release Date**: February 4, 2025

## 🎉 What's New

### ✨ Complete Workflow Documentation
- **Validated Quick Start Guide**: Step-by-step instructions tested end-to-end
- **Installation Scripts**: Enhanced macOS and Linux support with dependency management
- **Environment Setup**: Automated configuration of `EQDYNA2DCYCLEROOT` and PATH variables

### 📚 Comprehensive Documentation
- **README.md**: Complete rewrite with validated workflow commands
- **Academic Citations**: Properly formatted references to foundational research
- **Project Structure**: Clear organization and file descriptions
- **MIT License**: Open source licensing with separate LICENSE file

### 🔬 Validated Simulation Capabilities
Successfully tested complete workflow:
- **Multi-cycle Dynamics**: 10 earthquake cycles simulated
- **Realistic Results**: M5.6-M7.3+ earthquake magnitudes
- **Temporal Evolution**: Inter-event periods from 1-618 years
- **Automatic Visualization**: Post-processing with rupture dynamics plots

## 🛠️ Technical Improvements

### Enhanced Installation
```bash
# macOS with environment setup
./install.eqdyna.2dcycle.sh -e macos

# Ubuntu/Linux
./install.eqdyna.2dcycle.sh -m ubuntu
```

### Streamlined Workflow
```bash
# Complete workflow in 4 steps
export EQDYNA2DCYCLEROOT=$(pwd)
export PATH="$EQDYNA2DCYCLEROOT/bin:$EQDYNA2DCYCLEROOT/scripts:$PATH"
python3 scripts/create.newcase work/test_subei test.subei
cd work/test_subei && python3 case.setup && ./run.sh
```

### Robust Case Management
- Environment variable handling improvements
- Python path resolution fixes
- Automated file permission management
- Enhanced error handling

## 📊 Validation Results

### Test Case: Subei Fault System
- **Total Nodes**: 13,862
- **Total Elements**: 13,495
- **Model Range**: 136.86 × 106.03 km
- **Simulation Time**: Multi-cycle dynamics over 800+ years

### Representative Results
| Cycle | Magnitude | Max Slip (m) | Recurrence (years) |
|-------|-----------|--------------|-------------------|
| 1     | M7.33     | 3.95         | 618               |
| 2     | M5.69     | 0.26         | 141               |
| 3     | M5.79     | 0.31         | 28                |
| 4     | M5.61     | 0.22         | 1                 |
| 5     | M5.64     | 0.22         | 3                 |

## 📖 Academic References

### Core Methodology
> Duan, B., & Oglesby, D. D. (2006). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike‐slip faults. *Journal of Geophysical Research*, 111(B5), B05309. https://doi.org/10.1029/2005JB004138

### Applications
> Liu, D., Duan, B., Scharer, K., & Yule, D. (2022). Observation‐constrained multicycle dynamic models of the southern San Andreas and the northern San Jacinto faults. *Journal of Geophysical Research: Solid Earth*, 127(2), e2021JB023420. https://doi.org/10.1029/2021JB023420

## 🚀 Getting Started

1. **Clone the repository**
2. **Run installation**: `./install.eqdyna.2dcycle.sh -e macos`
3. **Create test case**: Follow the Quick Start guide in README.md
4. **View results**: Plots automatically generated in `aPlots/` directory

## 📧 Support

- **Website**: https://seismotamu.wixsite.com/emlam
- **Contact**: dliu@ig.utexas.edu, bduan@tamu.edu
- **Institution**: Earthquake Modeling Lab @ Texas A&M University

---

**Download**: Available on GitHub  
**License**: MIT License  
**Platform Support**: macOS, Linux, Unix-like systems