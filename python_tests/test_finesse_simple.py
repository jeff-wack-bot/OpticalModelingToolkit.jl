#!/usr/bin/env python3
"""
Simple test for Julia-Finesse interoperability
Usage: mamba run -n mtkoptical python test_finesse_simple.py
"""

import numpy as np
import sys
import unittest
import finesse
from juliacall import Main as jl


class TestJuliaFinesseSimple(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures."""
        try:
            jl.seval('import Pkg; Pkg.activate(".")')
        except:
            pass
    
    def test_finesse_import_and_version(self):
        """Test that Finesse can be imported and version accessed."""
        self.assertIsInstance(finesse.__version__, str)
        print(f"Finesse version: {finesse.__version__}")
        
        # Pass version to Julia
        jl.finesse_version = finesse.__version__
        version_from_julia = jl.seval("finesse_version")
        self.assertEqual(version_from_julia, finesse.__version__)
    
    def test_finesse_mirror_creation(self):
        """Test creating Finesse mirror and accessing from Julia."""
        from finesse.components.mirror import Mirror
        
        # Create mirror in Python
        mirror = Mirror(name="test_mirror", R=0.99, T=0.01)
        self.assertAlmostEqual(mirror.R.value, 0.99)
        self.assertAlmostEqual(mirror.T.value, 0.01)
        
        # Pass to Julia and access properties
        jl.mirror = mirror
        R_from_julia = jl.seval("mirror.R.value")
        T_from_julia = jl.seval("mirror.T.value")
        
        self.assertAlmostEqual(R_from_julia, 0.99)
        self.assertAlmostEqual(T_from_julia, 0.01)
        
        print(f"Mirror properties accessible from Julia: R={R_from_julia}, T={T_from_julia}")
    
    def test_optical_calculations_in_julia(self):
        """Test performing optical calculations in Julia with Python data."""
        # Create reflectivity data in Python
        R_values = np.array([0.9, 0.95, 0.99, 0.995, 0.999])
        
        # Pass to Julia
        jl.R_values = R_values
        
        # Perform cavity finesse calculation in Julia
        jl.seval("""
        using Statistics
        
        # Calculate cavity finesse: F = π√R / (1-R)
        finesse_values = [π * sqrt(R) / (1 - R) for R in R_values]
        
        # Calculate some statistics
        mean_finesse = mean(finesse_values)
        max_finesse = maximum(finesse_values)
        """)
        
        # Retrieve results
        finesse_results = np.array(jl.finesse_values)
        mean_finesse = jl.mean_finesse
        max_finesse = jl.max_finesse
        
        # Verify results
        self.assertEqual(len(finesse_results), len(R_values))
        self.assertTrue(np.all(finesse_results > 0))
        self.assertTrue(np.all(np.diff(finesse_results) > 0))  # Should be increasing
        self.assertGreater(max_finesse, mean_finesse)
        
        print(f"Cavity finesse calculations:")
        for R, F in zip(R_values, finesse_results):
            print(f"  R = {R:.3f} → Finesse = {F:.1f}")
        print(f"Mean finesse: {mean_finesse:.1f}")
        print(f"Max finesse: {max_finesse:.1f}")
    
    def test_data_array_exchange(self):
        """Test exchanging large arrays between Python and Julia."""
        # Create test data in Python
        frequencies = np.logspace(1, 5, 1000)  # 10 Hz to 100 kHz, 1000 points
        response = 1.0 / (1.0 + (frequencies / 1000.0)**2)  # Simple transfer function
        
        # Pass to Julia
        jl.frequencies = frequencies
        jl.response = response
        
        # Perform analysis in Julia
        jl.seval("""
        using Statistics
        
        # Find peak response and its frequency
        peak_idx = argmax(response)
        peak_freq = frequencies[peak_idx]
        peak_response = response[peak_idx]
        
        # Calculate 3dB bandwidth
        half_max = peak_response / sqrt(2)
        indices_3db = findall(r -> r >= half_max, response)
        if length(indices_3db) >= 2
            bandwidth_3db = frequencies[indices_3db[end]] - frequencies[indices_3db[1]]
        else
            bandwidth_3db = 0.0
        end
        
        # Calculate RMS
        rms_response = sqrt(mean(response.^2))
        """)
        
        # Retrieve results
        peak_freq = jl.peak_freq
        peak_response = jl.peak_response
        bandwidth_3db = jl.bandwidth_3db
        rms_response = jl.rms_response
        
        # Verify results make sense
        self.assertGreater(peak_response, 0.5)  # Should be near 1.0 at low frequencies
        self.assertGreater(bandwidth_3db, 0)
        self.assertGreater(rms_response, 0)
        self.assertLess(rms_response, peak_response)
        
        print(f"Transfer function analysis:")
        print(f"  Peak frequency: {peak_freq:.1f} Hz")
        print(f"  Peak response: {peak_response:.3f}")
        print(f"  3dB bandwidth: {bandwidth_3db:.1f} Hz")
        print(f"  RMS response: {rms_response:.3f}")
    
    def test_julia_mathematical_functions(self):
        """Test using Julia's mathematical functions on Finesse-related data."""
        # Cavity parameters
        mirror_R = [0.9, 0.95, 0.99, 0.995]
        cavity_length = 1.0  # meters
        
        jl.mirror_R = mirror_R
        jl.cavity_length = cavity_length
        
        # Calculate optical parameters in Julia
        jl.seval("""
        using Statistics, LinearAlgebra
        
        # Calculate free spectral range
        c = 299792458.0  # m/s
        fsr = c / (2 * cavity_length)
        
        # Calculate finesse for each mirror
        finesse_vals = [π * sqrt(R) / (1 - R) for R in mirror_R]
        
        # Calculate cavity bandwidth
        bandwidth_vals = fsr ./ finesse_vals
        
        # Statistical analysis
        mean_finesse = mean(finesse_vals)
        std_finesse = std(finesse_vals)
        """)
        
        # Get results
        fsr = jl.fsr
        finesse_vals = np.array(jl.finesse_vals)
        bandwidth_vals = np.array(jl.bandwidth_vals)
        mean_finesse = jl.mean_finesse
        std_finesse = jl.std_finesse
        
        # Verify physics
        self.assertAlmostEqual(fsr, 299792458.0 / 2.0, places=0)  # ~150 MHz
        self.assertEqual(len(finesse_vals), len(mirror_R))
        self.assertEqual(len(bandwidth_vals), len(mirror_R))
        self.assertTrue(np.all(finesse_vals > 0))
        self.assertTrue(np.all(bandwidth_vals > 0))
        
        # Higher reflectivity should give higher finesse and lower bandwidth
        self.assertTrue(np.all(np.diff(finesse_vals) > 0))  # Increasing finesse
        self.assertTrue(np.all(np.diff(bandwidth_vals) < 0))  # Decreasing bandwidth
        
        print(f"Cavity analysis results:")
        print(f"  Free spectral range: {fsr/1e6:.1f} MHz")
        print(f"  Mean finesse: {mean_finesse:.1f} ± {std_finesse:.1f}")
        
        for i, R in enumerate(mirror_R):
            print(f"  R = {R:.3f}: F = {finesse_vals[i]:.1f}, BW = {bandwidth_vals[i]/1e3:.1f} kHz")


def run_simple_tests():
    """Run simple interoperability tests."""
    print("="*60)
    print("SIMPLE JULIA-FINESSE INTEROPERABILITY TESTS")
    print("="*60)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJuliaFinesseSimple)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*60)
    if result.wasSuccessful():
        print("✓ All simple Julia-Finesse interoperability tests passed!")
    else:
        print(f"✗ {len(result.failures + result.errors)} tests failed!")
        return False
    print("="*60)
    return True


if __name__ == "__main__":
    success = run_simple_tests()
    sys.exit(0 if success else 1)
