# Basic Python-Julia interoperability tests for MTKOptical

# Configure PythonCall to use the mtkoptical environment
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/home/jeff/miniforge3/envs/mtkoptical/bin/python"

using PythonCall

@testset "Basic Python-Julia Interoperability" begin
    
    @testset "Python Environment Setup" begin
        # Test importing standard library
        sys = pyimport("sys")
        @test !isnothing(sys)
        @test pyhasattr(sys, "version")
        
        # Test Python executable path
        python_exe = PythonCall.python_executable_path()
        @test occursin("mtkoptical", python_exe)
    end
    
    @testset "NumPy Integration" begin
        np = pyimport("numpy")
        
        # Test array creation and operations
        x = np.array([1, 2, 3, 4, 5])
        @test pyconvert(Int, x.sum()) == 15
        @test pyconvert(Float64, x.mean()) ≈ 3.0
        
        # Test multidimensional arrays
        A = np.random.rand(3, 3)
        A_julia = pyconvert(Array, A)
        @test size(A_julia) == (3, 3)
        
        # Test mathematical operations
        B = np.eye(3)
        C = np.dot(A, B)
        @test pyconvert(Bool, np.allclose(A, C))
    end
    
    @testset "Data Type Conversion" begin
        # Test Julia → Python conversions
        julia_array = [1.0, 2.0, 3.0]
        py_array = PythonCall.Py(julia_array)
        @test pyconvert(Vector{Float64}, py_array) == julia_array
        
        # Test Python → Julia conversions
        np = pyimport("numpy")
        py_data = np.array([10, 20, 30])
        julia_data = pyconvert(Vector{Int}, py_data)
        @test julia_data == [10, 20, 30]
    end
    
    @testset "Mathematical Operations" begin
        np = pyimport("numpy")
        
        # Test linear algebra operations
        A = np.reshape(np.array([1.0, 2.0, 3.0, 4.0]), (2, 2))
        eigenvals, eigenvecs = np.linalg.eig(A)
        
        @test length(pyconvert(Vector, eigenvals)) == 2
        eigenvecs_julia = pyconvert(Array, eigenvecs)
        @test size(eigenvecs_julia) == (2, 2)
        
        # Verify eigenvalue computation worked without errors
        @test !isnothing(eigenvals)
        @test !isnothing(eigenvecs)
    end
    
    @testset "String and IO Operations" begin
        # Test string operations
        py_str = PythonCall.pystr("Hello from Julia!")
        @test pyconvert(String, py_str) == "Hello from Julia!"
        
        # Test file-like operations
        io = pyimport("io")
        string_buffer = io.StringIO("test data\nline 2")
        line1 = string_buffer.readline()
        @test pyconvert(String, line1) == "test data\n"
    end
    
    @testset "Error Handling" begin
        # Test Python exception handling
        @test_throws Exception pyimport("nonexistent_module")
        
        # Test calling nonexistent Python function
        math = pyimport("math")
        @test_throws Exception math.nonexistent_function()
    end
end