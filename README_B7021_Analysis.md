# B7021 Angular Contact Ball Bearing Analysis Program

## Overview

This Python program implements a comprehensive analysis of B7021 angular contact ball bearings based on Hertz contact theory. It calculates maximum contact stress and deformation for both inner and outer ring contacts under combined axial and radial loading conditions.

## Features

### Core Analysis Capabilities
- **Hertz Contact Theory Implementation**: Full mathematical model for elastic contact
- **Combined Loading Analysis**: Handles both axial and radial forces simultaneously
- **Load Distribution Calculation**: Determines load sharing among individual balls
- **Contact Stress Analysis**: Calculates maximum contact stress for inner and outer rings
- **Contact Deformation Analysis**: Computes elastic deformation at contact points
- **Engineering Assessment**: Provides safety factor analysis and design recommendations

### Advanced Features
- **Parametric Input Interface**: Easy modification of bearing parameters and load conditions
- **Sensitivity Analysis**: Study the effect of parameter variations on bearing performance
- **Multiple Configuration Support**: JSON-based configuration files for different scenarios
- **Comprehensive Visualization**: Detailed plots and charts for analysis results
- **Detailed Reporting**: Professional technical reports with mathematical derivations

## Bearing Specifications (B7021)

| Parameter | Value | Unit |
|-----------|-------|------|
| Number of balls (Z) | 19 | - |
| Ball diameter (D) | 19.04 | mm |
| Inner raceway diameter (d_i) | 105.0 | mm |
| Outer raceway diameter (d_o) | 160.0 | mm |
| Pitch diameter (d_m) | 132.5 | mm |
| Initial contact angle (α₀) | 15.0 | degrees |
| Inner curvature coefficient (f_i) | 0.515 | - |
| Outer curvature coefficient (f_o) | 0.520 | - |

## Load Conditions (Default)

| Load Type | Value | Unit |
|-----------|-------|------|
| Axial force (F_a) | 900 | N |
| Radial force (F_r) | 900 | N |

## Installation and Requirements

### Required Python Packages
```bash
pip install numpy scipy matplotlib pandas
```

### Python Version
- Python 3.7 or higher

## Usage

### Basic Analysis
```python
python3 b7021_bearing_analysis.py
```

This runs the complete analysis with default parameters and generates:
- Detailed console output with results
- Visualization plots (`b7021_bearing_analysis_results.png`)
- Technical report (`B7021_Bearing_Analysis_Report.txt`)

### Parametric Analysis
```python
python3 parametric_b7021_analysis.py
```

This performs sensitivity analysis and generates:
- Configuration files for different scenarios
- Sensitivity analysis plots
- Parametric analysis report

### Testing and Validation
```python
python3 test_b7021_analysis.py
```

Runs comprehensive validation tests to ensure calculation accuracy.

## File Structure

### Core Analysis Files
- `b7021_bearing_analysis.py` - Main analysis program
- `parametric_b7021_analysis.py` - Parametric interface and sensitivity analysis
- `test_b7021_analysis.py` - Validation tests

### Configuration Files (Generated)
- `b7021_base_config.json` - Default B7021 configuration
- `b7021_high_load_config.json` - High load scenario
- `b7021_light_load_config.json` - Light load scenario
- `b7021_ceramic_config.json` - Ceramic bearing material properties

### Output Files (Generated)
- `b7021_bearing_analysis_results.png` - Main analysis visualization
- `B7021_Bearing_Analysis_Report.txt` - Detailed technical report
- `sensitivity_*.png` - Sensitivity analysis plots
- `Parametric_Sensitivity_Report.txt` - Sensitivity analysis report

## Mathematical Model

### Hertz Contact Theory

The program implements the complete Hertz contact theory for elastic contact between curved surfaces:

#### Contact Deformation
```
δ = (9Q²/(16E*²(ρₓ + ρᵧ)))^(1/3)
```

#### Maximum Contact Stress
```
σₘₐₓ = 3Q/(2πab)
```

#### Contact Ellipse Semi-axes
```
a = (3QKₑ/(4E*(ρₓ + ρᵧ)))^(1/3)
b = (3QEₑ/(4E*(ρₓ + ρᵧ)))^(1/3)
```

Where:
- Q = Normal load on contact
- E* = Effective elastic modulus = E/(2(1-ν²))
- ρₓ, ρᵧ = Principal curvatures
- Kₑ, Eₑ = Complete elliptic integrals
- a, b = Semi-major and semi-minor axes of contact ellipse

### Curvature Calculations

#### Inner Ring Contact
```
ρₓᵢ = 1/r + 1/Rᵢ
ρᵧᵢ = 1/r
```

#### Outer Ring Contact
```
ρₓₒ = 1/r + 1/Rₒ
ρᵧₒ = 1/r
```

Where:
- r = Ball radius
- Rᵢ = Inner raceway radius = fᵢ × r
- Rₒ = Outer raceway radius = fₒ × r

## Sample Results

### Typical Analysis Output
```
BEARING PARAMETERS:
Number of balls (Z): 19
Ball diameter (D): 19.04 mm
Initial contact angle (α₀): 15.0°

LOAD CONDITIONS:
Axial force (F_a): 900 N
Radial force (F_r): 900 N
Total equivalent load: 3591.9 N
Maximum ball load: 283.2 N

CONTACT ANALYSIS RESULTS:
Inner Ring Contact:
  Maximum contact stress (σₙₐₓ): 3705.8 MPa
  Contact deformation (δ): 0.20 μm
  Contact area: 0.115 mm²

Outer Ring Contact:
  Maximum contact stress (σₙₐₓ): 3694.0 MPa
  Contact deformation (δ): 0.20 μm
  Contact area: 0.115 mm²

ENGINEERING ASSESSMENT:
Maximum contact stress in bearing: 3705.8 MPa
Safety factor: 0.5
⚠ Contact stress is high - consider design modifications
```

## Validation and Testing

The program includes comprehensive validation tests that verify:
- Parameter initialization correctness
- Mathematical model implementation
- Curvature calculations
- Contact deformation and stress calculations
- Load distribution analysis
- Geometric consistency
- Material property calculations

All tests pass validation, ensuring reliable calculations.

## Engineering Applications

### Design Analysis
- Bearing selection and verification
- Load capacity evaluation
- Safety factor assessment
- Performance optimization

### Research Applications
- Parameter sensitivity studies
- Material property effects
- Geometric optimization
- Load distribution analysis

### Quality Control
- Manufacturing tolerance analysis
- Performance validation
- Failure analysis support

## Limitations and Assumptions

### Analysis Limitations
- Static analysis (dynamic effects not considered)
- Ideal bearing geometry assumed
- Uniform material properties
- Lubrication effects not included
- Temperature effects not considered
- No consideration of manufacturing tolerances

### Model Assumptions
- Elastic contact behavior
- Hertzian contact conditions
- Linear load distribution approximation
- Perfect surface conditions

## Future Enhancements

### Potential Improvements
- Dynamic load analysis
- Thermal effects inclusion
- Lubrication modeling
- Fatigue life prediction
- Manufacturing tolerance effects
- Non-linear load distribution
- Multi-bearing system analysis

## References

1. Hertz, H. (1881). "On the contact of elastic solids"
2. Harris, T.A. & Kotzalas, M.N. (2007). "Rolling Bearing Analysis, 5th Edition"
3. ISO 281:2007 "Rolling bearings - Dynamic load ratings and rating life"
4. Palmgren, A. (1959). "Ball and Roller Bearing Engineering"

## Contact and Support

For technical questions or issues with the analysis program, please refer to the validation tests and mathematical documentation provided.

## License

This program is developed for educational and research purposes in bearing analysis and dynamics.

---

**Note**: This analysis program provides theoretical calculations based on Hertz contact theory. For critical applications, always validate results with experimental data and consult bearing manufacturers' specifications.