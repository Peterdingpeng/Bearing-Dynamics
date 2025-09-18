#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test module for B7021 Angular Contact Ball Bearing Analysis
Validates calculations and ensures reasonable results

Author: Bearing Dynamics Analysis
Last modified: 2024
"""

import numpy as np
import unittest
from b7021_bearing_analysis import (
    B7021BearingParameters, 
    HertzContactTheory, 
    LoadDistributionAnalysis, 
    BearingContactAnalysis
)

class TestB7021BearingAnalysis(unittest.TestCase):
    """Test class for B7021 bearing analysis"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.bearing_params = B7021BearingParameters()
        self.hertz_theory = HertzContactTheory(self.bearing_params)
        self.load_analysis = LoadDistributionAnalysis(self.bearing_params)
        self.contact_analysis = BearingContactAnalysis(self.bearing_params)
    
    def test_bearing_parameters(self):
        """Test bearing parameter initialization"""
        params = self.bearing_params
        
        # Check basic dimensions
        self.assertEqual(params.Z, 19)
        self.assertAlmostEqual(params.D * 1000, 19.04, places=2)  # Convert to mm
        self.assertAlmostEqual(params.d_i * 1000, 105.0, places=1)
        self.assertAlmostEqual(params.d_o * 1000, 160.0, places=1)
        self.assertAlmostEqual(params.alpha_0 * 180 / np.pi, 15.0, places=1)
        
        # Check calculated parameters
        self.assertAlmostEqual(params.r * 1000, 9.52, places=2)
        self.assertAlmostEqual(params.d_m * 1000, 132.5, places=1)
        
        # Check load conditions
        self.assertEqual(params.F_a, 900)
        self.assertEqual(params.F_r, 900)
    
    def test_curvature_calculations(self):
        """Test curvature parameter calculations"""
        rho_xi, rho_yi, rho_xo, rho_yo = self.hertz_theory.calculate_curvature_parameters()
        
        # All curvatures should be positive
        self.assertGreater(rho_xi, 0)
        self.assertGreater(rho_yi, 0)
        self.assertGreater(rho_xo, 0)
        self.assertGreater(rho_yo, 0)
        
        # Check reasonable magnitudes (order of 10² m⁻¹)
        self.assertGreater(rho_xi, 100)
        self.assertLess(rho_xi, 1000)
        self.assertGreater(rho_xo, 100)
        self.assertLess(rho_xo, 1000)
    
    def test_elliptical_parameters(self):
        """Test elliptical contact parameter calculations"""
        rho_x, rho_y = 300, 100  # Example curvature values
        F_rho, K_e, E_e = self.hertz_theory.calculate_elliptical_parameters(rho_x, rho_y)
        
        # F_rho should be between 0 and 1
        self.assertGreaterEqual(F_rho, 0)
        self.assertLessEqual(F_rho, 1)
        
        # Elliptic integrals should be positive
        self.assertGreater(K_e, 0)
        self.assertGreater(E_e, 0)
        
        # For this case, should be approximately π/2
        self.assertGreater(K_e, np.pi/4)
        self.assertLess(K_e, np.pi)
    
    def test_contact_deformation(self):
        """Test contact deformation calculations"""
        Q = 300  # Example load in N
        rho_x, rho_y = 300, 100  # Example curvatures
        
        delta, sigma_max, a, b = self.hertz_theory.calculate_contact_deformation(Q, rho_x, rho_y)
        
        # All results should be positive
        self.assertGreater(delta, 0)
        self.assertGreater(sigma_max, 0)
        self.assertGreater(a, 0)
        self.assertGreater(b, 0)
        
        # Check reasonable magnitudes
        self.assertLess(delta, 1e-5)  # Less than 10 μm
        self.assertGreater(sigma_max, 1e6)  # Greater than 1 MPa
        self.assertLess(sigma_max, 10e9)  # Less than 10 GPa
        self.assertLess(a, 1e-3)  # Less than 1 mm
        self.assertLess(b, 1e-3)  # Less than 1 mm
    
    def test_load_distribution(self):
        """Test load distribution calculations"""
        load_results = self.load_analysis.calculate_load_distribution()
        
        # Check that we have results for all balls
        self.assertEqual(len(load_results['Q_balls']), self.bearing_params.Z)
        self.assertEqual(len(load_results['ball_positions']), self.bearing_params.Z)
        
        # All loads should be non-negative
        self.assertTrue(np.all(load_results['Q_balls'] >= 0))
        
        # Total load should be reasonable
        self.assertGreater(load_results['Q_total'], 1000)  # Should be > 1000 N
        self.assertLess(load_results['Q_total'], 10000)   # Should be < 10000 N
        
        # Maximum ball load should be reasonable
        self.assertGreater(load_results['Q_max'], 100)    # Should be > 100 N
        self.assertLess(load_results['Q_max'], 1000)      # Should be < 1000 N
    
    def test_complete_analysis(self):
        """Test complete bearing analysis"""
        results = self.contact_analysis.analyze_contact_stress_and_deformation()
        
        # Check that all required keys are present
        required_keys = ['bearing_params', 'load_conditions', 'inner_ring_contact', 
                        'outer_ring_contact', 'load_distribution']
        for key in required_keys:
            self.assertIn(key, results)
        
        # Check inner ring results
        inner = results['inner_ring_contact']
        self.assertGreater(inner['sigma_max'], 0)
        self.assertGreater(inner['delta'], 0)
        self.assertGreater(inner['contact_area_a'], 0)
        self.assertGreater(inner['contact_area_b'], 0)
        
        # Check outer ring results
        outer = results['outer_ring_contact']
        self.assertGreater(outer['sigma_max'], 0)
        self.assertGreater(outer['delta'], 0)
        self.assertGreater(outer['contact_area_a'], 0)
        self.assertGreater(outer['contact_area_b'], 0)
        
        # Contact stress should be in reasonable range (1000-5000 MPa for heavily loaded bearings)
        self.assertGreater(inner['sigma_max'], 1000)
        self.assertLess(inner['sigma_max'], 10000)
        self.assertGreater(outer['sigma_max'], 1000)
        self.assertLess(outer['sigma_max'], 10000)
    
    def test_material_properties(self):
        """Test material property calculations"""
        # Check effective elastic modulus
        expected_E_eff = self.bearing_params.E / (2 * (1 - self.bearing_params.nu**2))
        self.assertAlmostEqual(self.hertz_theory.E_eff, expected_E_eff, places=0)
        
        # Check that it's reasonable for steel
        self.assertGreater(self.hertz_theory.E_eff, 100e9)  # > 100 GPa
        self.assertLess(self.hertz_theory.E_eff, 200e9)     # < 200 GPa
    
    def test_geometric_consistency(self):
        """Test geometric consistency of bearing parameters"""
        params = self.bearing_params
        
        # Outer diameter should be larger than inner diameter
        self.assertGreater(params.d_o, params.d_i)
        
        # Pitch diameter should be between inner and outer
        self.assertGreater(params.d_m, params.d_i)
        self.assertLess(params.d_m, params.d_o)
        
        # Ball diameter should be reasonable relative to raceway diameters
        self.assertLess(params.D, (params.d_o - params.d_i) / 2)
        
        # Curvature coefficients should be reasonable
        self.assertGreater(params.f_i, 0.5)
        self.assertLess(params.f_i, 1.0)
        self.assertGreater(params.f_o, 0.5)
        self.assertLess(params.f_o, 1.0)

def run_validation_tests():
    """Run validation tests and print results"""
    print("Running B7021 Bearing Analysis Validation Tests")
    print("=" * 55)
    
    # Create test suite
    test_suite = unittest.TestLoader().loadTestsFromTestCase(TestB7021BearingAnalysis)
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    # Print summary
    print("\nTest Summary:")
    print("-" * 20)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    if result.failures:
        print("\nFailures:")
        for test, traceback in result.failures:
            print(f"- {test}: {traceback}")
    
    if result.errors:
        print("\nErrors:")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback}")
    
    if len(result.failures) == 0 and len(result.errors) == 0:
        print("\n✓ All tests passed successfully!")
        return True
    else:
        print("\n✗ Some tests failed. Please review the implementation.")
        return False

def main():
    """Main testing function"""
    success = run_validation_tests()
    
    if success:
        print("\nValidation completed successfully!")
        print("The B7021 bearing analysis implementation is working correctly.")
    else:
        print("\nValidation failed!")
        print("Please review and fix the issues identified in the tests.")

if __name__ == "__main__":
    main()