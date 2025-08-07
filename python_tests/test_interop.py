#!/usr/bin/env python3
"""
Test script for Julia-Python interoperability
Usage: mamba run -n mtkoptical python test_interop.py
"""

import numpy as np
import sys
import unittest
from juliacall import Main as jl


class TestJuliaPythonInterop(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Activate the MTKOptical Julia project if it exists
        try:
            jl.seval('import Pkg; Pkg.activate(".")')
        except:
            pass  # Continue if project activation fails
    
    def test_basic_julia_operations(self):
        """Test basic Julia function calls and operations."""
        # Test basic arithmetic
        result = jl.sqrt(16)
        self.assertAlmostEqual(result, 4.0)
        
        # Test Julia version access
        version = jl.VERSION
        self.assertIsNotNone(version)
        print(f"Using Julia version: {version}")
    
    def test_julia_arrays_and_linear_algebra(self):
        """Test Julia array operations and linear algebra."""
        # Create Julia array and perform operations
        jl.seval("using LinearAlgebra")
        
        # Test array creation and basic operations
        jl.seval("test_array = [1, 2, 3, 4, 5]")
        julia_sum = jl.sum(jl.test_array)
        self.assertEqual(julia_sum, 15)
        
        # Test matrix operations
        jl.seval("A = rand(3, 3)")
        jl.seval("B = I(3)")  # Identity matrix
        jl.seval("C = A * B")
        
        # Verify A * I = A
        are_equal = jl.seval("A ≈ C")
        self.assertTrue(are_equal)
    
    def test_data_exchange(self):
        """Test bidirectional data exchange between Python and Julia."""
        # Python → Julia data transfer
        python_data = np.array([1.0, 2.0, 3.0, 4.0])
        jl.py_data = python_data
        
        # Perform Julia operations on Python data
        julia_result = jl.seval("sum(py_data)")
        expected_sum = np.sum(python_data)
        self.assertAlmostEqual(julia_result, expected_sum)
        
        # Julia → Python data transfer
        jl.seval("julia_data = collect(1:10)")
        python_result = np.asarray(jl.julia_data)
        expected = np.arange(1, 11)
        np.testing.assert_array_equal(python_result, expected)
    
    def test_statistical_operations(self):
        """Test statistical computations across languages."""
        # Generate data in Julia
        jl.seval("using Statistics")
        jl.seval("data = randn(1000)")
        
        # Get statistics from Julia
        julia_mean = jl.mean(jl.data)
        julia_std = jl.std(jl.data)
        
        # Convert to Python and verify
        python_data = np.asarray(jl.data)
        python_mean = np.mean(python_data)
        python_std = np.std(python_data, ddof=1)  # Julia uses n-1 denominator
        
        # Should be approximately equal (same data, same computations)
        self.assertAlmostEqual(julia_mean, python_mean, places=10)
        self.assertAlmostEqual(julia_std, python_std, places=10)
    
    def test_string_operations(self):
        """Test string handling between Julia and Python."""
        # Python string to Julia
        python_str = "Hello from Python!"
        jl.py_string = python_str
        
        # Julia string operations
        julia_result = jl.seval("uppercase(py_string)")
        self.assertEqual(julia_result, python_str.upper())
        
        # Julia string to Python
        jl.seval('julia_string = "Hello from Julia!"')
        python_result = str(jl.julia_string)
        self.assertEqual(python_result, "Hello from Julia!")
    
    def test_function_definition_and_calling(self):
        """Test defining functions in Julia and calling from Python."""
        # Define a Julia function
        jl.seval("""
        function quadratic_formula(a, b, c)
            discriminant = b^2 - 4*a*c
            if discriminant < 0
                return nothing
            end
            sqrt_disc = sqrt(discriminant)
            x1 = (-b + sqrt_disc) / (2*a)
            x2 = (-b - sqrt_disc) / (2*a)
            return (x1, x2)
        end
        """)
        
        # Call Julia function from Python
        # Solve x^2 - 5x + 6 = 0 (should give x = 2, 3)
        result = jl.quadratic_formula(1, -5, 6)
        x1, x2 = sorted([result[0], result[1]])  # Sort to ensure order
        
        self.assertAlmostEqual(x1, 2.0, places=10)
        self.assertAlmostEqual(x2, 3.0, places=10)
    
    def test_error_handling(self):
        """Test error handling across language boundaries."""
        # Test Julia error propagation to Python
        with self.assertRaises(Exception):
            jl.seval("error(\"This is a Julia error\")")
        
        # Test calling undefined Julia function
        with self.assertRaises(AttributeError):
            jl.nonexistent_function()
    
    def test_performance_critical_operation(self):
        """Test a performance-critical operation that benefits from Julia."""
        # Define a computational kernel in Julia
        jl.seval("""
        function matrix_multiply_sum(n)
            A = rand(n, n)
            B = rand(n, n)
            C = A * B
            return sum(C)
        end
        """)
        
        # Call from Python
        result = jl.matrix_multiply_sum(100)
        
        # Verify result is reasonable (sum of 100x100 random matrix products)
        self.assertIsInstance(result, (int, float))
        self.assertGreater(result, 0)  # Should be positive for random positive matrices


def run_comprehensive_tests():
    """Run all interoperability tests with detailed output."""
    print("="*60)
    print("PYTHON → JULIA INTEROPERABILITY TESTS")
    print("="*60)
    
    # Create test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(TestJuliaPythonInterop)
    
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*60)
    if result.wasSuccessful():
        print("✓ All Python → Julia interoperability tests passed!")
    else:
        print(f"✗ {len(result.failures + result.errors)} tests failed!")
        return False
    print("="*60)
    return True


if __name__ == "__main__":
    success = run_comprehensive_tests()
    sys.exit(0 if success else 1)
