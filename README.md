# OpticalModelingToolkit

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeff-wack-bot.github.io/OpticalModelingToolkit.jl/dev/)
[![Build Status](https://github.com/jeff-wack-bot/OpticalModelingToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeff-wack-bot/OpticalModelingToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

An optical domain for [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl).
It provides mirrors, beamsplitters, free space, lasers and photodetectors on
the plane-wave complex envelope. They compose with mechanics, electronics and
control in a single model, and are analyzed with the SciML stack: DC
operating points, linearized transfer functions, time-domain transients, and
exact frequency responses with every propagation delay kept transcendental.

```julia
eqs = [link(laser.p1, m1.p2); link(m1.p1, space.p1); link(space.p2, m2.p1); link(m2.p2, pd.p1)]
@named cavity = System(eqs, t; systems = [laser, m1, space, m2, pd])
```

See the [documentation](https://jeff-wack-bot.github.io/OpticalModelingToolkit.jl/dev/),
in particular the [design notes](https://jeff-wack-bot.github.io/OpticalModelingToolkit.jl/dev/design/).
