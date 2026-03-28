using CR3BPTools,
    NFCR3BP,
    StaticArrays,
    LinearAlgebra,
    DifferentialEquations,
    Roots,
    JLD2,
    NonlinearSolve,
    GLMakie,
    LaTeXStrings

include("../util/plotting/plotting.jl")

# ============================================================
# Three 3D trajectories side by side in one figure
# ============================================================

function plot_trio_3d(orbs::Vector{<:orbit_struct},
                      colors::Vector;
                      png_path = nothing)

    set_theme!(figure_settings(:small))
    fig = Figure(size=(1400, 480), figure_padding=(20, 20, 70, 30))

    az, el = -0.2π, 0.1π

    for (col, (orb, color)) in enumerate(zip(orbs, colors))
        traj = orb.syntraj
        x = getindex.(traj, 1)
        y = getindex.(traj, 2)
        z = getindex.(traj, 3)

        ax = Axis3(fig[1, col];
            xlabel       = L"x\ \mathrm{[N.D]}",
            ylabel       = L"y\ \mathrm{[N.D]}",
            zlabel       = L"z\ \mathrm{[N.D]}",
            title        = orb.name,
            aspect       = :data,
            azimuth      = az,
            elevation    = el,
            xlabeloffset = 40,
            ylabeloffset = 40,
            zlabeloffset = 80,
            zticks       = WilkinsonTicks(3),
        )

        lines!(ax, x, y, z; color=color, linewidth=1.5)
    end

    if !isnothing(png_path)
        save(png_path, fig; px_per_unit=2)
        println("Saved: $png_path")
    end

    return fig
end


## ── Load data ─────────────────────────────────────────────────────────────────
println("Loading data...")
L1_data_path = "data/Earth-Moon/Resonant/L1/Order16.jld2"

transformations = load(L1_data_path, "transformations")
hamiltonians    = load(L1_data_path, "hamiltonians")
derivatives     = load(L1_data_path, "derivatives")
eoms_AA         = derivatives.AA_flow

## ── Solve for halo actions ────────────────────────────────────────────────────
H = 0.95
θ2dot = eoms_AA[4]
θ3dot = eoms_AA[6]

function residual!(F, x, H_given)
    I2, I3 = x[1], x[2]
    state = SVector(0, 0, I2, pi/2, I3, 0)
    F[1] = evaluatepolynomial(θ2dot, state)
    F[2] = evaluatepolynomial(hamiltonians.action_angle, state) - H_given
end

sol = solve(NonlinearProblem(residual!, [0.4, 0.4], H),
            NewtonRaphson(); reltol=1e-12, abstol=1e-12)
I2_halo, I3_halo = sol.u

## ── Initial conditions ────────────────────────────────────────────────────────
# Lissajous:   θ2 = 0
# Quasihalo 1: θ2 = +π/2
# Quasihalo 2: θ2 = -π/2  (sign flip of 4th element)
tspan = LinRange(0, 100, 2000)

lissajous_AA  = SVector(0, 0, I2_halo + 0.05,       0,     I3_halo, 0)
quasihalo1_AA = SVector(0, 0, I2_halo + 0.1,  pi/2, I3_halo, 0)
quasihalo2_AA = SVector(0, 0, I2_halo + 0.1, -pi/2, I3_halo, 0)

## ── Propagate ─────────────────────────────────────────────────────────────────
lissajous_AAtraj  = propagate(eoms_polynomial, lissajous_AA,  tspan, eoms_AA)
println("Lissajous propagated")
quasihalo1_AAtraj = propagate(eoms_polynomial, quasihalo1_AA, tspan, eoms_AA)
println("Quasihalo 1 propagated")
quasihalo2_AAtraj = propagate(eoms_polynomial, quasihalo2_AA, tspan, eoms_AA)
println("Quasihalo 2 propagated")

## ── Build orbit structs ───────────────────────────────────────────────────────
lissajous_orbit = orbit_struct("Lissajous orbit",
    tspan,
    actionangle2synodic(lissajous_AAtraj,  transformations),
    actionangle2normalform(lissajous_AAtraj,  transformations),
    lissajous_AAtraj)

quasihalo1_orbit = orbit_struct("Northern quasihalo orbit",
    tspan,
    actionangle2synodic(quasihalo1_AAtraj, transformations),
    actionangle2normalform(quasihalo1_AAtraj, transformations),
    quasihalo1_AAtraj)

quasihalo2_orbit = orbit_struct("Southern quasihalo orbit",
    tspan,
    actionangle2synodic(quasihalo2_AAtraj, transformations),
    actionangle2normalform(quasihalo2_AAtraj, transformations),
    quasihalo2_AAtraj)

## ── Plot & save ───────────────────────────────────────────────────────────────
plot_trio_3d(
    [quasihalo1_orbit, lissajous_orbit, quasihalo2_orbit],
    [RGBf(0.25, 0.45, 0.75), RGBf(0.75, 0.28, 0.1), RGBf(0.8, 0.6, 0.1)];
    png_path = "script/nf_figs/trio_traj.png")
