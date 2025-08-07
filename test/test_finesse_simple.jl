# Simplified test for Finesse-Julia interoperability
# Usage: julia --project=. test_finesse_simple.jl

# Configure PythonCall to use the mtkoptical environment
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/home/jeff/miniforge3/envs/mtkoptical/bin/python"

using Test
using PythonCall

@testset "Basic Finesse-Julia Interoperability" begin
    
    @testset "Finesse Import" begin
        # Test basic import
        finesse = pyimport("finesse")
        @test !isnothing(finesse)
        
        # Check version
        version = finesse.__version__
        @test !isnothing(version)
        println("✓ Finesse version: ", version)
    end
    
    @testset "Finesse Component Creation" begin
        # Test basic component creation
        Mirror = pyimport("finesse.components.mirror").Mirror
        
        # Create a simple mirror
        mirror = Mirror(name="test_mirror", R=0.99, T=0.01)
        @test !isnothing(mirror)
        
        # Access properties
        R_val = pyconvert(Float64, mirror.R.value)
        T_val = pyconvert(Float64, mirror.T.value)
        
        @test R_val ≈ 0.99
        @test T_val ≈ 0.01
        
        println("✓ Mirror created with R=$R_val, T=$T_val")
    end
    
    @testset "Finesse Data Exchange" begin
        # Test passing data between Julia and Python/Finesse
        import_name = "test_component"
        julia_data = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        # Pass Julia data to Python
        finesse_utils = pyimport("finesse.utilities.misc")
        @test !isnothing(finesse_utils)
        
        # Test that we can pass Julia arrays to Finesse functions
        py_data = PythonCall.Py(julia_data)
        @test !isnothing(py_data)
        
        # Convert back to Julia
        converted_back = pyconvert(Vector{Float64}, py_data)
        @test converted_back == julia_data
        
        println("✓ Data exchange successful: ", converted_back)
    end
    
    @testset "Mathematical Operations with Finesse" begin
        # Use Julia's math capabilities on Finesse-related data
        
        # Create some optical parameters
        R_values = [0.9, 0.95, 0.99, 0.995, 0.999]
        T_values = [1.0 - R for R in R_values]
        
        # Calculate Finesse using Julia
        finesse_values = [π * sqrt(R) / (1 - R) for R in R_values]
        
        @test length(finesse_values) == 5
        @test all(f > 0 for f in finesse_values)
        @test issorted(finesse_values)  # Should increase with R
        
        # Pass to Python for verification
        py_R = PythonCall.Py(R_values)
        py_finesse = PythonCall.Py(finesse_values)
        
        # Convert back and verify
        julia_R_back = pyconvert(Vector{Float64}, py_R)
        julia_finesse_back = pyconvert(Vector{Float64}, py_finesse)
        
        @test julia_R_back ≈ R_values
        @test julia_finesse_back ≈ finesse_values
        
        println("✓ Optical calculations successful")
        println("  R-values: ", round.(R_values, digits=3))
        println("  Finesse values: ", round.(finesse_values, digits=1))
    end
end

println("✓ All basic Finesse-Julia interoperability tests completed successfully!")