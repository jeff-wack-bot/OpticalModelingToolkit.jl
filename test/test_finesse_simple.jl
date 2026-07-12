# # Finesse interoperability
#
# A smoke test for the Finesse–Julia bridge. [Finesse](https://finesse.ifosim.org/)
# is a Python interferometer simulator; OMTK reaches it through
# [PythonCall](https://juliapy.github.io/PythonCall.jl/). Finesse is installed by
# `CondaPkg.toml` at the project root (via pip), so this file needs no path
# configuration — it runs as part of `] test`.

using Test
using PythonCall

# ## Importing Finesse
#
# The Python package must import and report a version:

@testset "Finesse Import" begin
    finesse = pyimport("finesse")
    @test !isnothing(finesse)

    version = finesse.__version__
    @test !isnothing(version)
    println("✓ Finesse version: ", version)
end

# ## Creating a component
#
# A Finesse `Mirror` is created with reflectivity `R` and transmissivity `T`,
# and those values round-trip back through `pyconvert`:

@testset "Finesse Component Creation" begin
    Mirror = pyimport("finesse.components.mirror").Mirror

    R_val = 0.99
    T_val = 0.01

    mirror = Mirror(name = "test_mirror", R = R_val, T = T_val)
    @test !isnothing(mirror)

    ## properties are Python-side; convert them back to Float64 to compare
    @test pyconvert(Float64, mirror.R.value) ≈ R_val
    @test pyconvert(Float64, mirror.T.value) ≈ T_val

    println("✓ Mirror created with R=$R_val, T=$T_val")
end
