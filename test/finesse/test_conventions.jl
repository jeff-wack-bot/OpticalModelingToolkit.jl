# # Conventions: one surface at a time
#
# The mapping in `common.jl` (an OMT surface is a Finesse surface turned
# around, with `phi = -deg(ϕ)`; real transmission; √W amplitudes) is pinned here on single
# components, before anything is built on it. Each case drives one port and
# reads the complex amplitude leaving every port. OMT reads the exact
# scattering form (`to_scattering` at `f = 0`); Finesse reads `ad` detectors
# normalized to the field arriving at the driven port.

const RTϕ = [
    (0.36, 0.64, 0.0), (0.36, 0.64, deg2rad(10)), (0.99, 0.01, 0.3),
    (0.7, 0.2, -1.1), (1.0, 0.0, 0.2), (0.0, 1.0, 0.0),
]

# A one-port probe of `comp`: a laser into `comp.<drive>`, every other port
# closed with `Vacuum`. Returns `S[out]`, the amplitude leaving each port per
# unit amplitude entering `drive`.

function omt_column(comp, nports, drive)
    @named src = Laser()
    terms = [Symbol(:v, k) for k in 1:nports]
    vacs = [Vacuum(; name = terms[k]) for k in 1:nports if k != drive]
    vports = [k for k in 1:nports if k != drive]
    port(k) = getproperty(comp, Symbol(:p, k))
    eqs = vcat(link(src.p1, port(drive)), [link(port(k), v.p1) for (k, v) in zip(vports, vacs)]...)
    @named probe = System(eqs, t; systems = [src, comp, vacs...])
    parts = to_scattering(probe)
    G = omt_tf(parts, 0.0)
    iu = findfirst(==(:src), parts.u)
    names = Dict(k => (k == drive ? :src : terms[k]) for k in 1:nports)
    return [G[findfirst(==(names[k]), parts.y), iu] for k in 1:nports]
end

function finesse_column(katline, port, nports, drive)
    ads = join(("ad o$k $(port(k)).o 0" for k in 1:nports), "\n")
    sol = finesse_run(
        """
        laser l1 P=1
        $katline
        space s0 l1.p1 $(port(drive)) L=1
        ad in $(port(drive)).i 0
        $ads
        noxaxis()
        """
    )
    ain = reading(sol, "in")
    return [reading(sol, "o$k") / ain for k in 1:nports]
end

@testset "mirror scattering, both faces" begin
    for (R, T, ϕ) in RTϕ, drive in (1, 2)
        @named m = Mirror(; R, T, ϕ)
        omt = omt_column(m, 2, drive)
        fin = finesse_column(kat_mirror("m"; R, T, ϕ), k -> mport("m", k), 2, drive)
        @test omt ≈ fin rtol = 1.0e-12 atol = 1.0e-14
    end
end

@testset "beamsplitter scattering, all four ports" begin
    for (R, T, ϕ) in RTϕ, drive in 1:4
        @named b = Beamsplitter(; R, T, ϕ)
        omt = omt_column(b, 4, drive)
        fin = finesse_column(kat_bs("b"; R, T, ϕ), k -> bsport("b", k), 4, drive)
        @test omt ≈ fin rtol = 1.0e-12 atol = 1.0e-14
    end
end

# The laser's amplitude and phase, and the `ad` scale: a laser of power `P` and
# phase `ϕ` straight into a detector must read `√P e^{iϕ}` in both codes, and
# a `space` of any length must leave the carrier phase alone.

@testset "laser amplitude, phase, and detector units" begin
    for (P, ϕ, L) in ((1.0, 0.0, 1.0), (2.0, deg2rad(30), 1.234), (0.5, -2.0, 1.0e3))
        @named src = Laser(A = sqrt(P), ϕ = ϕ)
        @named sp = FreeSpace(; L)
        @named pd = Photodetector()
        @named line = System(vcat(link(src.p1, sp.p1), link(sp.p2, pd.p1)), t; systems = [src, sp, pd])
        sol = omt_dc(line)
        fin = finesse_run(
            """
            laser l1 P=$P phase=$(rad2deg(ϕ))
            mirror m R=0 T=1 L=0
            space s0 l1.p1 m.p1 L=$L
            ad a m.p1.i 0
            pd P m.p1.i
            noxaxis()
            """
        )
        @test cfield(sol, pd.p1.i) ≈ reading(fin, "a") rtol = 1.0e-12
        @test sol[pd.P] ≈ real(reading(fin, "P")) rtol = 1.0e-12
    end
end
