#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parametric Interface for B7021 Angular Contact Ball Bearing Analysis

This module provides a parametric input interface allowing users to:
1. Modify bearing parameters
2. Change load conditions
3. Adjust material properties
4. Run sensitivity analysis

Author: Bearing Dynamics Analysis
Last modified: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict, List, Tuple
import json
from b7021_bearing_analysis import (
    B7021BearingParameters, 
    HertzContactTheory, 
    LoadDistributionAnalysis, 
    BearingContactAnalysis
)

# Set matplotlib backend for headless environment
plt.switch_backend('Agg')

class ParametricBearingParameters(B7021BearingParameters):
    """Extended bearing parameters class with parametric input capability"""
    
    def __init__(self, config_dict: Dict = None):
        """Initialize with optional configuration dictionary"""
        super().__init__()
        
        if config_dict:
            self.load_from_dict(config_dict)
    
    def load_from_dict(self, config: Dict):
        """Load parameters from configuration dictionary"""
        
        # Bearing geometry
        if 'Z' in config:
            self.Z = config['Z']
        if 'D_mm' in config:
            self.D = config['D_mm'] * 1e-3  # Convert mm to m
            self.r = self.D / 2
        if 'd_i_mm' in config:
            self.d_i = config['d_i_mm'] * 1e-3  # Convert mm to m
        if 'd_o_mm' in config:
            self.d_o = config['d_o_mm'] * 1e-3  # Convert mm to m
        if 'alpha_0_deg' in config:
            self.alpha_0 = config['alpha_0_deg'] * np.pi / 180  # Convert deg to rad
        if 'f_i' in config:
            self.f_i = config['f_i']
        if 'f_o' in config:
            self.f_o = config['f_o']
            
        # Recalculate dependent parameters
        self.d_m = (self.d_i + self.d_o) / 2
        self.R_i = self.f_i * self.r
        self.R_o = self.f_o * self.r
        
        # Material properties
        if 'E_GPa' in config:
            self.E = config['E_GPa'] * 1e9  # Convert GPa to Pa
        if 'nu' in config:
            self.nu = config['nu']
        if 'rho_kg_m3' in config:
            self.rho = config['rho_kg_m3']
            
        # Load conditions
        if 'F_a' in config:
            self.F_a = config['F_a']
        if 'F_r' in config:
            self.F_r = config['F_r']
    
    def to_dict(self) -> Dict:
        """Export parameters to dictionary"""
        return {
            'Z': self.Z,
            'D_mm': self.D * 1000,
            'd_i_mm': self.d_i * 1000,
            'd_o_mm': self.d_o * 1000,
            'd_m_mm': self.d_m * 1000,
            'alpha_0_deg': self.alpha_0 * 180 / np.pi,
            'f_i': self.f_i,
            'f_o': self.f_o,
            'E_GPa': self.E / 1e9,
            'nu': self.nu,
            'rho_kg_m3': self.rho,
            'F_a': self.F_a,
            'F_r': self.F_r
        }
    
    def save_to_json(self, filename: str):
        """Save parameters to JSON file"""
        with open(filename, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
    
    def load_from_json(self, filename: str):
        """Load parameters from JSON file"""
        with open(filename, 'r') as f:
            config = json.load(f)
        self.load_from_dict(config)

class ParametricAnalysis:
    """Parametric analysis class for sensitivity studies"""
    
    def __init__(self, base_params: ParametricBearingParameters):
        self.base_params = base_params
    
    def sensitivity_analysis(self, parameter_name: str, variation_range: List[float]) -> Dict:
        """Perform sensitivity analysis on a single parameter"""
        
        results = {
            'parameter_values': [],
            'inner_stress': [],
            'outer_stress': [],
            'inner_deformation': [],
            'outer_deformation': [],
            'max_stress': [],
            'total_load': []
        }
        
        print(f"Performing sensitivity analysis on {parameter_name}")
        print(f"Variation range: {variation_range}")
        print("-" * 50)
        
        for i, value in enumerate(variation_range):
            # Create modified parameters
            config = self.base_params.to_dict()
            config[parameter_name] = value
            
            modified_params = ParametricBearingParameters(config)
            analyzer = BearingContactAnalysis(modified_params)
            
            # Run analysis (suppress output)
            import sys
            from io import StringIO
            
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            
            try:
                analysis_results = analyzer.analyze_contact_stress_and_deformation()
                
                # Store results
                results['parameter_values'].append(value)
                results['inner_stress'].append(analysis_results['inner_ring_contact']['sigma_max'])
                results['outer_stress'].append(analysis_results['outer_ring_contact']['sigma_max'])
                results['inner_deformation'].append(analysis_results['inner_ring_contact']['delta'])
                results['outer_deformation'].append(analysis_results['outer_ring_contact']['delta'])
                results['max_stress'].append(max(analysis_results['inner_ring_contact']['sigma_max'],
                                                analysis_results['outer_ring_contact']['sigma_max']))
                results['total_load'].append(analysis_results['load_conditions']['Q_total'])
                
            except Exception as e:
                print(f"Error at {parameter_name}={value}: {e}")
                # Use NaN for failed calculations
                results['parameter_values'].append(value)
                results['inner_stress'].append(np.nan)
                results['outer_stress'].append(np.nan)
                results['inner_deformation'].append(np.nan)
                results['outer_deformation'].append(np.nan)
                results['max_stress'].append(np.nan)
                results['total_load'].append(np.nan)
            
            finally:
                sys.stdout = old_stdout
            
            # Progress indicator
            progress = (i + 1) / len(variation_range) * 100
            print(f"Progress: {progress:.1f}% - {parameter_name}={value}")
        
        return results
    
    def multi_parameter_analysis(self, parameters: Dict[str, List[float]]) -> Dict:
        """Perform analysis with multiple parameter variations"""
        
        all_results = {}
        
        for param_name, param_range in parameters.items():
            print(f"\nAnalyzing parameter: {param_name}")
            results = self.sensitivity_analysis(param_name, param_range)
            all_results[param_name] = results
        
        return all_results
    
    def plot_sensitivity_results(self, results: Dict, parameter_name: str, save_path: str = None):
        """Plot sensitivity analysis results"""
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle(f'Sensitivity Analysis: {parameter_name}', fontsize=16, fontweight='bold')
        
        param_values = results['parameter_values']
        
        # Plot 1: Contact Stress vs Parameter
        axes[0, 0].plot(param_values, results['inner_stress'], 'o-', label='Inner Ring', linewidth=2, markersize=6)
        axes[0, 0].plot(param_values, results['outer_stress'], 's-', label='Outer Ring', linewidth=2, markersize=6)
        axes[0, 0].set_xlabel(parameter_name)
        axes[0, 0].set_ylabel('Maximum Contact Stress (MPa)')
        axes[0, 0].set_title('Contact Stress Sensitivity')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Contact Deformation vs Parameter
        axes[0, 1].plot(param_values, results['inner_deformation'], 'o-', label='Inner Ring', linewidth=2, markersize=6)
        axes[0, 1].plot(param_values, results['outer_deformation'], 's-', label='Outer Ring', linewidth=2, markersize=6)
        axes[0, 1].set_xlabel(parameter_name)
        axes[0, 1].set_ylabel('Contact Deformation (μm)')
        axes[0, 1].set_title('Contact Deformation Sensitivity')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Maximum Stress vs Parameter
        axes[1, 0].plot(param_values, results['max_stress'], 'ro-', linewidth=2, markersize=6)
        axes[1, 0].set_xlabel(parameter_name)
        axes[1, 0].set_ylabel('Maximum Stress in Bearing (MPa)')
        axes[1, 0].set_title('Overall Maximum Stress')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Add safety limit line
        axes[1, 0].axhline(y=2000, color='orange', linestyle='--', linewidth=2, label='Typical Yield (2000 MPa)')
        axes[1, 0].legend()
        
        # Plot 4: Total Load vs Parameter
        axes[1, 1].plot(param_values, results['total_load'], 'go-', linewidth=2, markersize=6)
        axes[1, 1].set_xlabel(parameter_name)
        axes[1, 1].set_ylabel('Total Equivalent Load (N)')
        axes[1, 1].set_title('Load Sensitivity')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Sensitivity plot saved to: {save_path}")
        else:
            plt.savefig(f'/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/sensitivity_{parameter_name}.png', 
                       dpi=300, bbox_inches='tight')
            print(f"Sensitivity plot saved to: sensitivity_{parameter_name}.png")
        
        plt.close()

def create_example_configurations():
    """Create example configuration files"""
    
    # Base B7021 configuration
    base_config = {
        "Z": 19,
        "D_mm": 19.04,
        "d_i_mm": 105.0,
        "d_o_mm": 160.0,
        "alpha_0_deg": 15.0,
        "f_i": 0.515,
        "f_o": 0.520,
        "E_GPa": 210.0,
        "nu": 0.3,
        "rho_kg_m3": 7850,
        "F_a": 900,
        "F_r": 900
    }
    
    # High load configuration
    high_load_config = base_config.copy()
    high_load_config.update({
        "F_a": 1500,
        "F_r": 1500
    })
    
    # Light load configuration
    light_load_config = base_config.copy()
    light_load_config.update({
        "F_a": 500,
        "F_r": 500
    })
    
    # Different material configuration (ceramic bearing)
    ceramic_config = base_config.copy()
    ceramic_config.update({
        "E_GPa": 320.0,    # Silicon nitride
        "nu": 0.27,
        "rho_kg_m3": 3200
    })
    
    # Save configurations
    configs = {
        'b7021_base_config.json': base_config,
        'b7021_high_load_config.json': high_load_config,
        'b7021_light_load_config.json': light_load_config,
        'b7021_ceramic_config.json': ceramic_config
    }
    
    for filename, config in configs.items():
        filepath = f'/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/{filename}'
        with open(filepath, 'w') as f:
            json.dump(config, f, indent=4)
        print(f"Created configuration file: {filename}")

def run_parametric_analysis_examples():
    """Run example parametric analyses"""
    
    print("Running Parametric Analysis Examples")
    print("=" * 40)
    
    # Create base parameters
    base_params = ParametricBearingParameters()
    analyzer = ParametricAnalysis(base_params)
    
    # Example 1: Load sensitivity analysis
    print("\n1. Load Sensitivity Analysis")
    print("-" * 30)
    
    load_range = np.linspace(500, 1500, 11)  # Vary axial load from 500 to 1500 N
    load_results = analyzer.sensitivity_analysis('F_a', load_range.tolist())
    analyzer.plot_sensitivity_results(load_results, 'Axial Load F_a (N)')
    
    # Example 2: Ball diameter sensitivity
    print("\n2. Ball Diameter Sensitivity Analysis")
    print("-" * 40)
    
    diameter_range = np.linspace(15, 23, 9)  # Vary ball diameter from 15 to 23 mm
    diameter_results = analyzer.sensitivity_analysis('D_mm', diameter_range.tolist())
    analyzer.plot_sensitivity_results(diameter_results, 'Ball Diameter D (mm)')
    
    # Example 3: Contact angle sensitivity
    print("\n3. Contact Angle Sensitivity Analysis")
    print("-" * 40)
    
    angle_range = np.linspace(10, 25, 8)  # Vary contact angle from 10 to 25 degrees
    angle_results = analyzer.sensitivity_analysis('alpha_0_deg', angle_range.tolist())
    analyzer.plot_sensitivity_results(angle_results, 'Contact Angle α₀ (degrees)')
    
    # Create summary report
    create_sensitivity_report(load_results, diameter_results, angle_results)

def create_sensitivity_report(load_results: Dict, diameter_results: Dict, angle_results: Dict):
    """Create a comprehensive sensitivity analysis report"""
    
    report_filename = '/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/Parametric_Sensitivity_Report.txt'
    
    with open(report_filename, 'w', encoding='utf-8') as f:
        f.write("B7021 BEARING PARAMETRIC SENSITIVITY ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. ANALYSIS OVERVIEW\n")
        f.write("-" * 20 + "\n")
        f.write("This report presents the results of parametric sensitivity analysis\n")
        f.write("performed on the B7021 angular contact ball bearing.\n\n")
        f.write("Parameters analyzed:\n")
        f.write("• Axial load (F_a): 500 - 1500 N\n")
        f.write("• Ball diameter (D): 15 - 23 mm\n")
        f.write("• Contact angle (α₀): 10 - 25 degrees\n\n")
        
        f.write("2. LOAD SENSITIVITY RESULTS\n")
        f.write("-" * 30 + "\n")
        
        # Find min/max stress values for load variation
        min_stress = min(load_results['max_stress'])
        max_stress = max(load_results['max_stress'])
        load_sensitivity = (max_stress - min_stress) / min_stress * 100
        
        f.write(f"Load range analyzed: {min(load_results['parameter_values'])} - {max(load_results['parameter_values'])} N\n")
        f.write(f"Minimum maximum stress: {min_stress:.1f} MPa\n")
        f.write(f"Maximum maximum stress: {max_stress:.1f} MPa\n")
        f.write(f"Stress sensitivity: {load_sensitivity:.1f}% variation\n\n")
        
        f.write("3. BALL DIAMETER SENSITIVITY RESULTS\n")
        f.write("-" * 40 + "\n")
        
        # Find min/max stress values for diameter variation
        min_stress = min(diameter_results['max_stress'])
        max_stress = max(diameter_results['max_stress'])
        diameter_sensitivity = (max_stress - min_stress) / min_stress * 100
        
        f.write(f"Diameter range analyzed: {min(diameter_results['parameter_values'])} - {max(diameter_results['parameter_values'])} mm\n")
        f.write(f"Minimum maximum stress: {min_stress:.1f} MPa\n")
        f.write(f"Maximum maximum stress: {max_stress:.1f} MPa\n")
        f.write(f"Stress sensitivity: {diameter_sensitivity:.1f}% variation\n\n")
        
        f.write("4. CONTACT ANGLE SENSITIVITY RESULTS\n")
        f.write("-" * 40 + "\n")
        
        # Find min/max stress values for angle variation
        min_stress = min(angle_results['max_stress'])
        max_stress = max(angle_results['max_stress'])
        angle_sensitivity = (max_stress - min_stress) / min_stress * 100
        
        f.write(f"Angle range analyzed: {min(angle_results['parameter_values'])} - {max(angle_results['parameter_values'])} degrees\n")
        f.write(f"Minimum maximum stress: {min_stress:.1f} MPa\n")
        f.write(f"Maximum maximum stress: {max_stress:.1f} MPa\n")
        f.write(f"Stress sensitivity: {angle_sensitivity:.1f}% variation\n\n")
        
        f.write("5. SENSITIVITY RANKING\n")
        f.write("-" * 25 + "\n")
        
        sensitivities = [
            ('Axial Load', load_sensitivity),
            ('Ball Diameter', diameter_sensitivity),
            ('Contact Angle', angle_sensitivity)
        ]
        
        # Sort by sensitivity (highest first)
        sensitivities.sort(key=lambda x: x[1], reverse=True)
        
        f.write("Parameters ranked by sensitivity (highest to lowest):\n")
        for i, (param, sensitivity) in enumerate(sensitivities, 1):
            f.write(f"{i}. {param}: {sensitivity:.1f}% variation\n")
        
        f.write("\n6. DESIGN RECOMMENDATIONS\n")
        f.write("-" * 30 + "\n")
        f.write("Based on the sensitivity analysis:\n\n")
        
        if sensitivities[0][1] > 50:
            f.write(f"• {sensitivities[0][0]} shows high sensitivity ({sensitivities[0][1]:.1f}%)\n")
            f.write("  → This parameter requires careful control during design and operation\n")
        
        if sensitivities[-1][1] < 20:
            f.write(f"• {sensitivities[-1][0]} shows low sensitivity ({sensitivities[-1][1]:.1f}%)\n")
            f.write("  → This parameter has less impact on bearing performance\n")
        
        f.write("\n• Monitor parameters with high stress sensitivity closely\n")
        f.write("• Consider tolerance requirements based on sensitivity rankings\n")
        f.write("• Implement design margins for high-sensitivity parameters\n\n")
        
        f.write("Report generated by Parametric B7021 Bearing Analysis Program\n")
    
    print(f"Sensitivity analysis report saved to: {report_filename}")

def main():
    """Main parametric analysis function"""
    
    print("B7021 Parametric Bearing Analysis Interface")
    print("=" * 45)
    
    # Create example configurations
    print("\nCreating example configuration files...")
    create_example_configurations()
    
    # Run parametric analysis examples
    print("\nRunning parametric analysis examples...")
    run_parametric_analysis_examples()
    
    print("\nParametric analysis completed!")
    print("\nFiles generated:")
    print("- Configuration files: b7021_*_config.json")
    print("- Sensitivity plots: sensitivity_*.png")
    print("- Sensitivity report: Parametric_Sensitivity_Report.txt")
    
    # Demonstrate loading from configuration file
    print("\nDemonstrating configuration file usage...")
    demo_config_usage()

def demo_config_usage():
    """Demonstrate how to use configuration files"""
    
    print("\nConfiguration File Usage Example:")
    print("-" * 40)
    
    # Load high load configuration
    high_load_params = ParametricBearingParameters()
    high_load_params.load_from_json('/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/b7021_high_load_config.json')
    
    # Run analysis
    analyzer = BearingContactAnalysis(high_load_params)
    results = analyzer.analyze_contact_stress_and_deformation()
    
    print("High Load Configuration Results:")
    print(f"Axial Load: {results['load_conditions']['F_a']} N")
    print(f"Radial Load: {results['load_conditions']['F_r']} N")
    print(f"Max Inner Ring Stress: {results['inner_ring_contact']['sigma_max']:.1f} MPa")
    print(f"Max Outer Ring Stress: {results['outer_ring_contact']['sigma_max']:.1f} MPa")

if __name__ == "__main__":
    main()