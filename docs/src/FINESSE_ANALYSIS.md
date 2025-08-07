# Finesse 3 Architecture Analysis

Based on analysis of the Finesse 3 codebase structure and key components for developing MTKOptics.jl.

## Core Structure

### Model-Centric Design
The `Model` class (`src/finesse/model.py`) is the central container that manages all components, connections, and simulations. This uses NetworkX for graph-based component relationships and provides:
- Component registration and management
- Connection validation and network building
- Simulation orchestration
- Parameter management across the system

### Component Hierarchy
All optical elements inherit from well-defined base classes:
- **`Connector`** - Base for components with ports and connections
- **`Surface`** - Base for optical surfaces (mirrors, beamsplitters) 
- **`FrequencyGenerator`** - Base for sources (lasers, modulators)
- **`ModelElement`** - Root base class for all model elements

## Key Components Found

### Optical Elements
- **Mirror** (`mirror.py`) - Reflectivity, transmission, loss, curvature, tuning
- **Laser** (`laser.py`) - Power, frequency, phase control
- **Space** (`space.py`) - Propagation distances with refractive index
- **Beamsplitter** (`beamsplitter.py`) - Amplitude/phase splitting
- **Lens** (`lens.py`) - Focal length transformations
- **Modulator** (`modulator.py`) - Phase/amplitude modulation

### Detection and Measurement
- **PowerDetector** - Optical power measurements
- **Camera** - Spatial beam profiling
- **Cavity detectors** - Resonance measurements
- **Quantum noise detectors** - Shot noise, radiation pressure

### Mechanical Systems
- **Free masses** - Mechanical resonators
- **Suspension systems** - Pendulum dynamics
- **Radiation pressure** - Optical force coupling

### Electronics
- **Amplifiers** - Signal conditioning
- **Filters** - Transfer function processing
- **Signal processing** - Feedback and control

## Critical Design Patterns

### Port-Node System
Each component has ports (`p1`, `p2`, etc.) with directional nodes for beam paths:
```python
self._add_port("p1", NodeType.OPTICAL)
self.p1._add_node("i", NodeDirection.INPUT)
self.p1._add_node("o", NodeDirection.OUTPUT)
```
This enables:
- Automatic connection validation
- Bidirectional beam propagation
- Clear component interface definition

### Parameter System
Heavy use of decorators for parameter validation and units:
```python
@float_parameter("R", "Reflectivity", validate="_check_R", post_validate="check_rtl")
@float_parameter("phi", "Microscopic tuning", units="degrees")
@bool_parameter("misaligned", "Misaligns mirror reflection")
```
Features:
- Automatic validation and type checking
- Unit handling and conversion
- Parameter dependencies and constraints

### Cython Integration
Performance-critical code written in Cython (`.pyx` files):
- Matrix operations (`cymath/` directory)
- Simulation solvers (`simulations/*.pyx`)
- Component modal calculations (`components/modal/*.pyx`)
- Provides C-level performance for numerical kernels

### Symbolic Mathematics
Custom symbolic system for preprocessing equations:
- Symbolic parameter handling
- Automatic differentiation preparation
- Expression optimization before numerical evaluation

## Simulation Architecture

### Graph-Based Modeling
- Converts component network to mathematical graph representation
- Uses NetworkX for graph algorithms and analysis
- Enables automatic equation generation from network topology

### Sparse Matrix Methods
- Uses KLU sparse solver for large linear systems
- Optimized for interferometer network sparsity patterns
- Handles systems with thousands of optical modes

### Frequency Domain Calculations
- All calculations performed in frequency domain
- Multiple sideband support for modulated signals
- Hermite-Gaussian spatial mode expansion
- Complex amplitude representation throughout

### Solver Integration
Multiple solver backends:
- **BaseSolver** - Core solver interface
- **Sparse solvers** - KLU integration for large systems
- **HOM solvers** - Higher-order mode calculations
- **Specialized solvers** - Quantum noise, radiation pressure

## Implementation Insights for MTKOptics.jl

### Patterns to Adopt
1. **Port-node connection system** - Clean component interfaces
2. **Parameter decoration pattern** - Validation and units handling  
3. **Graph-based model representation** - Network topology management
4. **Separation of definition from simulation** - Components define structure, solvers handle numerics

### Patterns to Adapt
1. **Replace Cython with Julia native performance** - Leverage Julia's speed without compilation complexity
2. **Integrate MTK's symbolic preprocessing** - Use ModelingToolkit's symbolic capabilities
3. **Use SciML sparse solvers** - Leverage Julia's advanced solver ecosystem
4. **Adopt Julia's type system** - More sophisticated type handling than Python decorators

### Key Architectural Insight
Finesse separates component *definition* from *simulation*:
- Components define connectivity, parameters, and physical properties
- Separate simulation classes handle numerical solving and matrix operations  
- This clean separation works well with ModelingToolkit's approach

The codebase demonstrates mature patterns for optical simulation that can be effectively adapted to Julia's ecosystem while leveraging MTK's symbolic-numeric strengths and the broader SciML solver capabilities.

## Directory Structure Summary

```
src/finesse/
├── components/          # Optical/mechanical/electronic components
│   ├── modal/          # Cython implementations for performance
│   ├── mirror.py       # Mirror component definition
│   ├── laser.py        # Laser source definition
│   └── space.py        # Propagation space definition
├── simulations/        # Numerical solvers and simulation engines
│   ├── sparse/         # Sparse matrix solver integration
│   └── graph/          # Graph-based simulation tools
├── detectors/          # Measurement and detection components
├── cymath/             # Cython mathematical utilities
├── tracing/            # Beam propagation and ABCD matrix tools
└── model.py            # Central model container and orchestration
```

This analysis provides a solid foundation for implementing MTKOptics.jl with similar architectural principles adapted to Julia's strengths.