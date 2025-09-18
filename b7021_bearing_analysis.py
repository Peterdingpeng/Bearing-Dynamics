#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B7021 Angular Contact Ball Bearing Analysis Program
Based on Hertz Contact Theory

This program analyzes B7021 angular contact ball bearing under combined axial and radial loads,
calculating maximum contact stress and deformation using Hertz contact theory.

Author: Bearing Dynamics Analysis
Last modified: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, minimize_scalar
import pandas as pd
from typing import Dict, Tuple, List
import warnings

# Set matplotlib backend for headless environment
plt.switch_backend('Agg')

class B7021BearingParameters:
    """B7021 Angular Contact Ball Bearing Parameters"""
    
    def __init__(self):
        # Basic bearing geometry (from problem statement)
        self.Z = 19                    # Number of balls
        self.alpha_0 = 15 * np.pi / 180  # Initial contact angle (rad)
        self.D = 19.04e-3              # Ball diameter (m)
        self.r = self.D / 2            # Ball radius (m)
        self.d_i = 105e-3              # Inner raceway diameter (m)
        self.d_o = 160e-3              # Outer raceway diameter (m)
        self.d_m = (self.d_i + self.d_o) / 2  # Pitch diameter (m)
        self.f_i = 0.515               # Inner raceway curvature coefficient
        self.f_o = 0.52                # Outer raceway curvature coefficient
        
        # Calculated parameters
        self.R_i = self.f_i * self.r   # Inner raceway radius
        self.R_o = self.f_o * self.r   # Outer raceway radius
        
        # Material properties (typical bearing steel)
        self.E = 210e9                 # Young's modulus (Pa)
        self.nu = 0.3                  # Poisson's ratio
        self.rho = 7850                # Density (kg/m³)
        
        # Load conditions (from problem statement)
        self.F_a = 900                 # Axial force (N)
        self.F_r = 900                 # Radial force (N)

class HertzContactTheory:
    """Implementation of Hertz Contact Theory for Ball Bearings"""
    
    def __init__(self, bearing_params: B7021BearingParameters):
        self.params = bearing_params
        
        # Calculate effective elastic modulus
        self.E_eff = self.params.E / (2 * (1 - self.params.nu**2))
        
    def calculate_curvature_parameters(self) -> Tuple[float, float, float, float]:
        """Calculate curvature parameters for inner and outer contacts"""
        
        # Inner ring contact (corrected for proper curvature calculation)
        # For inner ring: ball curvature + raceway curvature (both positive)
        rho_xi = 1 / self.params.r + 1 / self.params.R_i
        rho_yi = 1 / self.params.r
        
        # Outer ring contact  
        # For outer ring: ball curvature + raceway curvature
        rho_xo = 1 / self.params.r + 1 / self.params.R_o
        rho_yo = 1 / self.params.r
        
        return rho_xi, rho_yi, rho_xo, rho_yo
    
    def calculate_elliptical_parameters(self, rho_x: float, rho_y: float) -> Tuple[float, float, float]:
        """Calculate elliptical contact parameters"""
        
        # Ensure positive curvatures
        rho_x = abs(rho_x)
        rho_y = abs(rho_y)
        
        # Curvature difference and sum
        F_rho = abs(rho_x - rho_y) / (rho_x + rho_y)
        
        # Elliptic integrals approximation
        if abs(F_rho) < 1e-6:  # Circular contact
            K_e = np.pi / 2
            E_e = np.pi / 2
        else:
            # First-order approximation for elliptic integrals
            K_e = np.pi / 2 * (1 + F_rho**2 / 4)
            E_e = np.pi / 2 * (1 - F_rho**2 / 4)
            
        return F_rho, K_e, E_e
    
    def calculate_contact_deformation(self, Q: float, rho_x: float, rho_y: float) -> Tuple[float, float, float]:
        """Calculate contact deformation and stress using Hertz theory"""
        
        # Ensure positive curvatures
        rho_x = abs(rho_x)
        rho_y = abs(rho_y)
        
        F_rho, K_e, E_e = self.calculate_elliptical_parameters(rho_x, rho_y)
        
        # Contact deformation
        delta = (9 * Q**2 / (16 * self.E_eff**2 * (rho_x + rho_y)))**(1/3)
        
        # Semi-axes of contact ellipse
        a = (3 * Q * K_e / (4 * self.E_eff * (rho_x + rho_y)))**(1/3)
        b = (3 * Q * E_e / (4 * self.E_eff * (rho_x + rho_y)))**(1/3)
        
        # Maximum contact stress (at center of contact)
        sigma_max = 3 * Q / (2 * np.pi * a * b)
        
        return delta, sigma_max, a, b

class LoadDistributionAnalysis:
    """Analysis of load distribution in angular contact ball bearing"""
    
    def __init__(self, bearing_params: B7021BearingParameters):
        self.params = bearing_params
        self.hertz = HertzContactTheory(bearing_params)
        
    def calculate_equilibrium_contact_angle(self, Q_total: float) -> float:
        """Calculate equilibrium contact angle under load"""
        
        def equilibrium_equation(alpha):
            # Simplified approach - contact angle increases with load
            # This is an approximation; full solution requires iterative analysis
            return alpha - self.params.alpha_0 - Q_total * 1e-6
            
        try:
            alpha_eq = fsolve(equilibrium_equation, self.params.alpha_0)[0]
            return max(self.params.alpha_0, min(alpha_eq, np.pi/3))  # Limit reasonable range
        except:
            return self.params.alpha_0
    
    def calculate_load_distribution(self) -> Dict:
        """Calculate load distribution among balls considering combined loading"""
        
        # Total equivalent load (simplified approach for combined radial and axial)
        Q_total = np.sqrt(self.params.F_r**2 + (self.params.F_a / np.sin(self.params.alpha_0))**2)
        
        # Load per ball (assuming equal distribution as first approximation)
        Q_ball_avg = Q_total / self.params.Z
        
        # Calculate load variation due to radial load
        ball_positions = np.linspace(0, 2*np.pi, self.params.Z, endpoint=False)
        
        # Load distribution considering radial load direction
        Q_balls = np.zeros(self.params.Z)
        for i, theta in enumerate(ball_positions):
            # Simplified load distribution - maximum at bottom (theta=3π/2)
            load_factor = 1 + 0.5 * np.cos(theta - 3*np.pi/2)  # Bottom ball carries more load
            Q_balls[i] = Q_ball_avg * load_factor
            
        # Ensure non-negative loads
        Q_balls = np.maximum(Q_balls, 0)
        
        # Normalize to maintain equilibrium
        Q_balls = Q_balls * Q_total / np.sum(Q_balls)
        
        return {
            'Q_total': Q_total,
            'Q_balls': Q_balls,
            'ball_positions': ball_positions,
            'max_loaded_ball': np.argmax(Q_balls),
            'Q_max': np.max(Q_balls)
        }

class BearingContactAnalysis:
    """Main class for bearing contact analysis"""
    
    def __init__(self, bearing_params: B7021BearingParameters):
        self.params = bearing_params
        self.hertz = HertzContactTheory(bearing_params)
        self.load_dist = LoadDistributionAnalysis(bearing_params)
        
    def analyze_contact_stress_and_deformation(self) -> Dict:
        """Perform complete contact analysis"""
        
        print("Starting B7021 Angular Contact Ball Bearing Analysis...")
        print("=" * 60)
        
        # Calculate load distribution
        load_results = self.load_dist.calculate_load_distribution()
        Q_max = load_results['Q_max']
        
        # Get curvature parameters
        rho_xi, rho_yi, rho_xo, rho_yo = self.hertz.calculate_curvature_parameters()
        
        print(f"Curvature Parameters:")
        print(f"Inner ring: ρ_xi = {rho_xi:.2e} m⁻¹, ρ_yi = {rho_yi:.2e} m⁻¹")
        print(f"Outer ring: ρ_xo = {rho_xo:.2e} m⁻¹, ρ_yo = {rho_yo:.2e} m⁻¹")
        print()
        
        # Analyze inner ring contact
        delta_i, sigma_max_i, a_i, b_i = self.hertz.calculate_contact_deformation(Q_max, rho_xi, rho_yi)
        
        # Analyze outer ring contact
        delta_o, sigma_max_o, a_o, b_o = self.hertz.calculate_contact_deformation(Q_max, rho_xo, rho_yo)
        
        results = {
            'bearing_params': {
                'Z': self.params.Z,
                'D': self.params.D * 1000,  # Convert to mm for display
                'd_i': self.params.d_i * 1000,
                'd_o': self.params.d_o * 1000,
                'd_m': self.params.d_m * 1000,
                'alpha_0_deg': self.params.alpha_0 * 180 / np.pi,
                'f_i': self.params.f_i,
                'f_o': self.params.f_o
            },
            'load_conditions': {
                'F_a': self.params.F_a,
                'F_r': self.params.F_r,
                'Q_total': load_results['Q_total'],
                'Q_max': Q_max
            },
            'inner_ring_contact': {
                'sigma_max': sigma_max_i / 1e6,  # Convert to MPa
                'delta': delta_i * 1e6,  # Convert to μm
                'contact_area_a': a_i * 1e3,  # Convert to mm
                'contact_area_b': b_i * 1e3,  # Convert to mm
                'contact_area': np.pi * a_i * b_i * 1e6  # Convert to mm²
            },
            'outer_ring_contact': {
                'sigma_max': sigma_max_o / 1e6,  # Convert to MPa
                'delta': delta_o * 1e6,  # Convert to μm
                'contact_area_a': a_o * 1e3,  # Convert to mm
                'contact_area_b': b_o * 1e3,  # Convert to mm
                'contact_area': np.pi * a_o * b_o * 1e6  # Convert to mm²
            },
            'load_distribution': load_results
        }
        
        return results
    
    def print_results(self, results: Dict):
        """Print detailed analysis results"""
        
        print("BEARING PARAMETERS:")
        print("-" * 30)
        params = results['bearing_params']
        print(f"Number of balls (Z): {params['Z']}")
        print(f"Ball diameter (D): {params['D']:.2f} mm")
        print(f"Inner raceway diameter (d_i): {params['d_i']:.1f} mm")
        print(f"Outer raceway diameter (d_o): {params['d_o']:.1f} mm")
        print(f"Pitch diameter (d_m): {params['d_m']:.1f} mm")
        print(f"Initial contact angle (α₀): {params['alpha_0_deg']:.1f}°")
        print(f"Inner curvature coefficient (f_i): {params['f_i']:.3f}")
        print(f"Outer curvature coefficient (f_o): {params['f_o']:.3f}")
        print()
        
        print("LOAD CONDITIONS:")
        print("-" * 30)
        loads = results['load_conditions']
        print(f"Axial force (F_a): {loads['F_a']:.0f} N")
        print(f"Radial force (F_r): {loads['F_r']:.0f} N")
        print(f"Total equivalent load: {loads['Q_total']:.1f} N")
        print(f"Maximum ball load: {loads['Q_max']:.1f} N")
        print()
        
        print("CONTACT ANALYSIS RESULTS:")
        print("-" * 30)
        
        inner = results['inner_ring_contact']
        outer = results['outer_ring_contact']
        
        print("Inner Ring Contact:")
        print(f"  Maximum contact stress (σ_max): {inner['sigma_max']:.1f} MPa")
        print(f"  Contact deformation (δ): {inner['delta']:.2f} μm")
        print(f"  Contact ellipse semi-axes: a = {inner['contact_area_a']:.3f} mm, b = {inner['contact_area_b']:.3f} mm")
        print(f"  Contact area: {inner['contact_area']:.3f} mm²")
        print()
        
        print("Outer Ring Contact:")
        print(f"  Maximum contact stress (σ_max): {outer['sigma_max']:.1f} MPa")
        print(f"  Contact deformation (δ): {outer['delta']:.2f} μm")
        print(f"  Contact ellipse semi-axes: a = {outer['contact_area_a']:.3f} mm, b = {outer['contact_area_b']:.3f} mm")
        print(f"  Contact area: {outer['contact_area']:.3f} mm²")
        print()
        
        # Engineering assessment
        print("ENGINEERING ASSESSMENT:")
        print("-" * 30)
        max_stress = max(inner['sigma_max'], outer['sigma_max'])
        typical_yield = 2000  # MPa for bearing steel
        safety_factor = typical_yield / max_stress
        
        print(f"Maximum contact stress in bearing: {max_stress:.1f} MPa")
        print(f"Typical bearing steel yield strength: {typical_yield} MPa")
        print(f"Safety factor: {safety_factor:.1f}")
        
        if safety_factor > 3:
            print("✓ Contact stress is within safe operating limits")
        elif safety_factor > 1.5:
            print("⚠ Contact stress is moderate - monitor operation")
        else:
            print("⚠ Contact stress is high - consider design modifications")
    
    def create_visualizations(self, results: Dict):
        """Create visualization plots"""
        
        # Set up the plotting style
        plt.style.use('default')
        fig = plt.figure(figsize=(15, 10))
        
        # Plot 1: Load distribution
        ax1 = plt.subplot(2, 3, 1)
        load_dist = results['load_distribution']
        ball_numbers = np.arange(1, self.params.Z + 1)
        plt.bar(ball_numbers, load_dist['Q_balls'], alpha=0.7, color='skyblue', edgecolor='black')
        plt.xlabel('Ball Number')
        plt.ylabel('Load (N)')
        plt.title('Load Distribution Among Balls')
        plt.grid(True, alpha=0.3)
        
        # Plot 2: Contact stress comparison
        ax2 = plt.subplot(2, 3, 2)
        stress_types = ['Inner Ring', 'Outer Ring']
        stress_values = [results['inner_ring_contact']['sigma_max'], 
                        results['outer_ring_contact']['sigma_max']]
        colors = ['lightcoral', 'lightgreen']
        bars = plt.bar(stress_types, stress_values, color=colors, alpha=0.7, edgecolor='black')
        plt.ylabel('Maximum Contact Stress (MPa)')
        plt.title('Contact Stress Comparison')
        plt.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, stress_values):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                    f'{value:.1f}', ha='center', va='bottom', fontweight='bold')
        
        # Plot 3: Contact deformation comparison
        ax3 = plt.subplot(2, 3, 3)
        deform_values = [results['inner_ring_contact']['delta'], 
                        results['outer_ring_contact']['delta']]
        bars = plt.bar(stress_types, deform_values, color=['orange', 'purple'], alpha=0.7, edgecolor='black')
        plt.ylabel('Contact Deformation (μm)')
        plt.title('Contact Deformation Comparison')
        plt.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, deform_values):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{value:.2f}', ha='center', va='bottom', fontweight='bold')
        
        # Plot 4: Contact area visualization
        ax4 = plt.subplot(2, 3, 4)
        contact_areas = [results['inner_ring_contact']['contact_area'], 
                        results['outer_ring_contact']['contact_area']]
        bars = plt.bar(stress_types, contact_areas, color=['gold', 'lightblue'], alpha=0.7, edgecolor='black')
        plt.ylabel('Contact Area (mm²)')
        plt.title('Contact Area Comparison')
        plt.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, contact_areas):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # Plot 5: Load distribution polar plot
        ax5 = plt.subplot(2, 3, 5, projection='polar')
        theta = load_dist['ball_positions']
        loads = load_dist['Q_balls']
        ax5.plot(theta, loads, 'o-', linewidth=2, markersize=6, color='red')
        ax5.fill(theta, loads, alpha=0.3, color='red')
        ax5.set_title('Polar Load Distribution', pad=20)
        ax5.set_ylim(0, max(loads) * 1.1)
        
        # Plot 6: Summary data table
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('tight')
        ax6.axis('off')
        
        summary_data = [
            ['Parameter', 'Inner Ring', 'Outer Ring'],
            ['Max Stress (MPa)', f"{results['inner_ring_contact']['sigma_max']:.1f}", 
             f"{results['outer_ring_contact']['sigma_max']:.1f}"],
            ['Deformation (μm)', f"{results['inner_ring_contact']['delta']:.2f}", 
             f"{results['outer_ring_contact']['delta']:.2f}"],
            ['Contact Area (mm²)', f"{results['inner_ring_contact']['contact_area']:.3f}", 
             f"{results['outer_ring_contact']['contact_area']:.3f}"],
            ['Semi-axis a (mm)', f"{results['inner_ring_contact']['contact_area_a']:.3f}", 
             f"{results['outer_ring_contact']['contact_area_a']:.3f}"],
            ['Semi-axis b (mm)', f"{results['inner_ring_contact']['contact_area_b']:.3f}", 
             f"{results['outer_ring_contact']['contact_area_b']:.3f}"]
        ]
        
        table = ax6.table(cellText=summary_data[1:], colLabels=summary_data[0], 
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.5)
        ax6.set_title('Contact Analysis Summary')
        
        plt.tight_layout()
        plt.savefig('/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/b7021_bearing_analysis_results.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualization saved as 'b7021_bearing_analysis_results.png'")
    
    def generate_report(self, results: Dict):
        """Generate detailed calculation report"""
        
        report_filename = '/home/runner/work/Bearing-Dynamics/Bearing-Dynamics/B7021_Bearing_Analysis_Report.txt'
        
        with open(report_filename, 'w', encoding='utf-8') as f:
            f.write("B7021 ANGULAR CONTACT BALL BEARING ANALYSIS REPORT\n")
            f.write("=" * 60 + "\n")
            f.write("Based on Hertz Contact Theory\n\n")
            
            f.write("1. BEARING SPECIFICATIONS\n")
            f.write("-" * 30 + "\n")
            params = results['bearing_params']
            f.write(f"Bearing Type: B7021 Angular Contact Ball Bearing\n")
            f.write(f"Number of balls (Z): {params['Z']}\n")
            f.write(f"Ball diameter (D): {params['D']:.2f} mm\n")
            f.write(f"Ball radius (r): {params['D']/2:.2f} mm\n")
            f.write(f"Inner raceway diameter (d_i): {params['d_i']:.1f} mm\n")
            f.write(f"Outer raceway diameter (d_o): {params['d_o']:.1f} mm\n")
            f.write(f"Pitch diameter (d_m): {params['d_m']:.1f} mm\n")
            f.write(f"Initial contact angle (α₀): {params['alpha_0_deg']:.1f}°\n")
            f.write(f"Inner raceway curvature coefficient (f_i): {params['f_i']:.3f}\n")
            f.write(f"Outer raceway curvature coefficient (f_o): {params['f_o']:.3f}\n\n")
            
            f.write("2. LOAD CONDITIONS\n")
            f.write("-" * 30 + "\n")
            loads = results['load_conditions']
            f.write(f"Axial force (F_a): {loads['F_a']:.0f} N\n")
            f.write(f"Radial force (F_r): {loads['F_r']:.0f} N\n")
            f.write(f"Total equivalent load: {loads['Q_total']:.1f} N\n")
            f.write(f"Maximum ball load: {loads['Q_max']:.1f} N\n\n")
            
            f.write("3. HERTZ CONTACT THEORY CALCULATIONS\n")
            f.write("-" * 40 + "\n")
            f.write("3.1 Mathematical Model:\n")
            f.write("The Hertz contact theory for elastic contact between two bodies:\n\n")
            f.write("Contact deformation: δ = (9Q²/(16E*²(ρₓ + ρᵧ)))^(1/3)\n")
            f.write("Maximum contact stress: σₘₐₓ = 3Q/(2πab)\n")
            f.write("Semi-axes: a = (3QKₑ/(4E*(ρₓ + ρᵧ)))^(1/3)\n")
            f.write("           b = (3QEₑ/(4E*(ρₓ + ρᵧ)))^(1/3)\n\n")
            f.write("Where:\n")
            f.write("Q = Normal load on contact\n")
            f.write("E* = Effective elastic modulus = E/(2(1-ν²))\n")
            f.write("ρₓ, ρᵧ = Principal curvatures\n")
            f.write("Kₑ, Eₑ = Complete elliptic integrals\n\n")
            
            f.write("3.2 Material Properties:\n")
            f.write(f"Young's modulus (E): {self.params.E/1e9:.0f} GPa\n")
            f.write(f"Poisson's ratio (ν): {self.params.nu:.1f}\n")
            f.write(f"Effective modulus (E*): {self.hertz.E_eff/1e9:.1f} GPa\n\n")
            
            f.write("4. CONTACT ANALYSIS RESULTS\n")
            f.write("-" * 35 + "\n")
            
            inner = results['inner_ring_contact']
            outer = results['outer_ring_contact']
            
            f.write("4.1 Inner Ring Contact:\n")
            f.write(f"Maximum contact stress (σₘₐₓ): {inner['sigma_max']:.1f} MPa\n")
            f.write(f"Contact deformation (δ): {inner['delta']:.2f} μm\n")
            f.write(f"Contact ellipse semi-major axis (a): {inner['contact_area_a']:.3f} mm\n")
            f.write(f"Contact ellipse semi-minor axis (b): {inner['contact_area_b']:.3f} mm\n")
            f.write(f"Contact area: {inner['contact_area']:.3f} mm²\n\n")
            
            f.write("4.2 Outer Ring Contact:\n")
            f.write(f"Maximum contact stress (σₘₐₓ): {outer['sigma_max']:.1f} MPa\n")
            f.write(f"Contact deformation (δ): {outer['delta']:.2f} μm\n")
            f.write(f"Contact ellipse semi-major axis (a): {outer['contact_area_a']:.3f} mm\n")
            f.write(f"Contact ellipse semi-minor axis (b): {outer['contact_area_b']:.3f} mm\n")
            f.write(f"Contact area: {outer['contact_area']:.3f} mm²\n\n")
            
            f.write("5. ENGINEERING EVALUATION\n")
            f.write("-" * 30 + "\n")
            max_stress = max(inner['sigma_max'], outer['sigma_max'])
            typical_yield = 2000  # MPa for bearing steel
            safety_factor = typical_yield / max_stress
            
            f.write(f"Maximum contact stress in bearing: {max_stress:.1f} MPa\n")
            f.write(f"Location of maximum stress: {'Inner ring' if inner['sigma_max'] > outer['sigma_max'] else 'Outer ring'}\n")
            f.write(f"Typical bearing steel yield strength: {typical_yield} MPa\n")
            f.write(f"Safety factor: {safety_factor:.1f}\n\n")
            
            if safety_factor > 3:
                f.write("Assessment: Contact stress is within safe operating limits.\n")
            elif safety_factor > 1.5:
                f.write("Assessment: Contact stress is moderate - monitor operation closely.\n")
            else:
                f.write("Assessment: Contact stress is high - consider design modifications.\n")
            
            f.write("\n6. RECOMMENDATIONS\n")
            f.write("-" * 20 + "\n")
            f.write("• Regular monitoring of bearing temperature and vibration\n")
            f.write("• Proper lubrication maintenance\n")
            f.write("• Consider load reduction if operating near stress limits\n")
            f.write("• Periodic inspection for wear and fatigue\n\n")
            
            f.write("7. ANALYSIS LIMITATIONS\n")
            f.write("-" * 25 + "\n")
            f.write("• Static analysis - dynamic effects not considered\n")
            f.write("• Assumes ideal bearing geometry\n")
            f.write("• Material properties assumed uniform\n")
            f.write("• Lubrication effects not included\n")
            f.write("• Temperature effects not considered\n\n")
            
            f.write("Report generated by B7021 Bearing Analysis Program\n")
            f.write("Based on Hertz Contact Theory\n")
        
        print(f"Detailed report saved as '{report_filename}'")

def main():
    """Main analysis function"""
    
    print("B7021 Angular Contact Ball Bearing Analysis")
    print("Based on Hertz Contact Theory")
    print("=" * 50)
    print()
    
    # Initialize bearing parameters
    bearing = B7021BearingParameters()
    
    # Create analysis object
    analyzer = BearingContactAnalysis(bearing)
    
    # Perform analysis
    results = analyzer.analyze_contact_stress_and_deformation()
    
    # Print results
    analyzer.print_results(results)
    
    # Create visualizations
    analyzer.create_visualizations(results)
    
    # Generate detailed report
    analyzer.generate_report(results)
    
    print("\nAnalysis completed successfully!")
    print("Files generated:")
    print("- b7021_bearing_analysis_results.png (visualization)")
    print("- B7021_Bearing_Analysis_Report.txt (detailed report)")

if __name__ == "__main__":
    main()