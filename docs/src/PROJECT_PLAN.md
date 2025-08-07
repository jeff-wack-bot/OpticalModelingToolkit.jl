# MTKOptics.jl Project Plan

## Vision
Create a Julia optical interferometer modeling package combining Finesse 3's proven optical concepts with ModelingToolkit.jl's symbolic-numeric framework for high-performance simulations.

## Core Strategy
- ✅ **Finesse 3 Analysis**: Cloned repository is located in "./finesse3". Currently this is not installed as a dev dependency, so it is just to look but not touch.
- ✅ **Python-Julia Interoperability**: Full bidirectional communication
- ✅ **Finesse Integration**: Can call Finesse from Julia for ground truth validation
- ✅ **Test Infrastructure**: Include comparisons against finesse

### Validation-Driven Development
**Use Finesse 3 as Ground Truth**: MTKOptical implementations will be validated against Finesse results
- Finesse provides proven optical simulation algorithms for comparison
- Julia can call Finesse directly to generate reference solutions
- Test-driven development ensures compatibility and correctness
- Focus on Julia→Python/Finesse calls for validation (primary workflow)

### Architecture Approach
**Hybrid Design**: Leverage ModelingToolkit's strengths while adding specialized optical mathematics
- The key goal here is to allow for simulation of quantum-electro-opto-mechanical systems.
- How can we make a system which is able to simulate hybrid quantum-classical systems as described in Wiseman and Milburn?
- Use MTK for: symbolic preprocessing, automatic differentiation, component composition, sparse matrix optimization
- Custom additions: optical field representations, Hermite-Gaussian modes, frequency-domain equations

## Key Technical Decisions

### Numerical Strategy
- **Matrix Operations**: Julia's sparse structures for interferometer networks
- **Solvers**: NonlinearSolve.jl with MTK symbolic acceleration for steady-state solutions
- **Performance**: Symbolic sparsity detection, type stability, SIMD vectorization

### Component Interface
```julia
@component function Mirror(; name, R=0.99, T=0.01, phi=0.0)
    @parameters R T phi
    ports = [OpticalPort(; name=:p1), OpticalPort(; name=:p2)]
    # Transfer matrix symbolic relationships
end
```

### Syntax Philosophy
**Native Julia** approach (not Finesse DSL clone):
```julia
sys = @optical_system begin
    L0 = Laser(P=1.0, λ=1064e-9)
    m1 = Mirror(R=0.99, T=0.01)
    cavity = Space(L=1.0)
    connect(L0.out, m1.p1)
end
```

**Question** Can the ModelingToolkit syntax be used exactly in order to maximize interoperability with that ecosystem?

## Implementation Priorities

### Phase 0: Infrastructure ✅ COMPLETED 
- ✅ **Environment Setup**: Python-Julia interoperability fully configured with miniforge3/mamba
- ✅ **Finesse Integration**: Finesse 3.0a33 installed and accessible from Julia (PyPI version)
- ✅ **Validation Framework**: Comprehensive test suites (30/30 Julia tests, 13/13 Python tests passing)
- ✅ **Test Organization**: Proper Julia package structure with working/failing test separation
- ✅ **Documentation**: Complete setup guide with troubleshooting and best practices
- ✅ **Ground Truth Capability**: Verified Julia→Finesse calls for validation workflows

### Phase 1: Foundation with Validation (CURRENT FOCUS)
- **Core Components**: Laser, Mirror, Space, BeamSplitter with MTK integration
- **Finesse Comparison Tests**: Each component validated against Finesse equivalent
- **OpticalField Type**: Complex amplitude representation compatible with Finesse
- **Connection System**: Port-based connections similar to Finesse architecture
- **Validation Workflow**: Automated tests comparing MTKOptical vs Finesse results

### Phase 2: Core Physics with Benchmarking
- **Hermite-Gaussian Modes**: Implementation validated against Finesse HOM calculations
- **Frequency-Domain Propagation**: Transfer matrix formalism with Finesse cross-checks
- **Steady-State Solving**: NonlinearSolve.jl integration with Finesse result validation
- **Performance Benchmarks**: Speed comparisons between MTKOptical and Finesse

### Phase 3: Advanced Features and Optimization
- **Nonlinear Crystals**: Advanced optical elements with validation
- **Quantum Noise Integration**: Quantum optics features beyond Finesse capabilities
- **Performance Optimization**: GPU acceleration while maintaining validation
- **Migration Tools**: Utilities for converting Finesse models to MTKOptical


## Development Workflow

### Primary Validation Pattern
```julia
# Infrastructure already established - example workflow:

# 1. Implement MTKOptical component
mtk_result = simulate(mtk_model)

# 2. Create equivalent Finesse model via Julia→Python (WORKING)
# Configure environment
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/home/username/miniforge3/envs/mtkoptical/bin/python"
using PythonCall

# Import Finesse and create reference
finesse = pyimport("finesse")
finesse_model = create_finesse_reference(mtk_model)
finesse_result = run_finesse_simulation(finesse_model)

# 3. Validate results
@test isapprox(mtk_result, finesse_result, rtol=1e-3)
```

### Test Infrastructure
- **Julia Package Tests**: `test/runtests.jl` with basic interop and Finesse integration
- **Python Test Suite**: `python_tests/` with bidirectional functionality verification  
- **Failing Test Analysis**: `failing_tests/` with post-mortem analysis for future debugging
- **All Tests Passing**: 30/30 Julia tests, 13/13 Python tests confirmed working

### Ground Truth Philosophy
- **Finesse as Oracle**: Use Finesse's proven algorithms as reference standard ✅ **IMPLEMENTED**
- **Iterative Refinement**: Improve MTKOptical until it matches Finesse results
- **Performance Focus**: Once validated, optimize for Julia's performance advantages
- **Extension Beyond**: Add new capabilities not available in Finesse

### Established Capabilities
- **Bidirectional Interoperability**: Zero-copy array transfers, efficient data exchange
- **Mathematical Integration**: Julia's advanced math on Finesse simulation data
- **Scientific Computing**: Optical physics calculations (cavity finesse, FSR, bandwidth)
- **Robust Error Handling**: Exception propagation across language boundaries
- **Performance Verified**: Large dataset processing (1000+ element arrays) working efficiently

This foundation ensures that MTKOptical development can proceed with confidence, maintaining the scientific rigor of Finesse while unlocking Julia's performance and ecosystem advantages.
