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
# Static plotting functions (PNG output, no animation)
# ============================================================

function plot_orbit_phasespace_static(orb::orbit_struct,
                                      transformations;
                                      png_path = nothing,
                                      color    = :royalblue)

    traj = orb.syntraj
    NF   = orb.NFtraj

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)
    z  = getindex.(traj, 3)

    q2 = getindex.(NF, 3)
    p2 = getindex.(NF, 4)
    q3 = getindex.(NF, 5)
    p3 = getindex.(NF, 6)

    xL1 = transformations.recentering[1]

    set_theme!(figure_settings(:large))
    fig = Figure(size=(1600, 800), figure_padding=(80, 20, 20, 20))

    ax3d = Axis3(fig[1:2, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]", zlabel=L"z [N.D]",
        aspect=:data,
        title=orb.name,
        azimuth   = π/3,
        elevation = 0.15π,
    )
    ax_q2p2 = Axis(fig[1, 2];
        xlabel=L"\bar{q}_2", ylabel=L"\bar{p}_2",
        title=L"\bar{q}_2\ \mathrm{vs}\ \bar{p}_2",
        aspect=DataAspect(),
    )
    ax_q3p3 = Axis(fig[2, 2];
        xlabel=L"\bar{q}_3", ylabel=L"\bar{p}_3",
        title=L"\bar{q}_3\ \mathrm{vs}\ \bar{p}_3",
        aspect=DataAspect(),
    )

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    lines!(ax3d,    x,  y,  z;  color=color, linewidth=2)
    lines!(ax_q2p2, q2, p2;     color=color, linewidth=1.5)
    lines!(ax_q3p3, q3, p3;     color=color, linewidth=1.5)

    scatter!(ax3d, [xL1], [0.0], [0.0]; color=:red, markersize=12)

    if !isnothing(png_path)
        save(png_path, fig; px_per_unit=2)
        println("Saved: $png_path")
    end

    return fig
end


function plot_orbit_phasespace_2d_static(orb::orbit_struct,
                                          transformations;
                                          png_path = nothing,
                                          color    = :royalblue)

    traj = orb.syntraj
    NF   = orb.NFtraj

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)

    q2 = getindex.(NF, 3)
    p2 = getindex.(NF, 4)
    q3 = getindex.(NF, 5)
    p3 = getindex.(NF, 6)

    xL1 = transformations.recentering[1]

    set_theme!(figure_settings(:large))
    fig = Figure(size=(1600, 800), figure_padding=(40, 20, 20, 20))

    ax_xy = Axis(fig[1:2, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]",
        title=orb.name,
        aspect=DataAspect(),
    )
    ax_q2p2 = Axis(fig[1, 2];
        xlabel=L"\bar{q}_2", ylabel=L"\bar{p}_2",
        title=L"\bar{q}_2\ \mathrm{vs}\ \bar{p}_2",
        aspect=DataAspect(),
    )
    ax_q3p3 = Axis(fig[2, 2];
        xlabel=L"\bar{q}_3", ylabel=L"\bar{p}_3",
        title=L"\bar{q}_3\ \mathrm{vs}\ \bar{p}_3",
        aspect=DataAspect(),
    )

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    lines!(ax_xy,   x,  y;   color=color, linewidth=2)
    lines!(ax_q2p2, q2, p2;  color=color, linewidth=1.5)
    lines!(ax_q3p3, q3, p3;  color=color, linewidth=1.5)

    scatter!(ax_xy, [xL1], [0.0]; color=:red, markersize=12)

    if !isnothing(png_path)
        save(png_path, fig; px_per_unit=2)
        println("Saved: $png_path")
    end

    return fig
end


function plot_orbit_actionangle_2d_static(orb::orbit_struct,
                                           transformations;
                                           png_path = nothing,
                                           color    = :royalblue)

    traj = orb.syntraj
    AA   = orb.AAtraj
    ts   = collect(orb.tspan)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)

    I2 = getindex.(AA, 3)
    θ2 = mod.(getindex.(AA, 4), 2π)
    I3 = getindex.(AA, 5)
    θ3 = mod.(getindex.(AA, 6), 2π)

    xL1 = transformations.recentering[1]

    set_theme!(figure_settings(:large))
    fig = Figure(size=(1600, 900), figure_padding=(40, 20, 20, 20))

    ax_xy = Axis(fig[1:4, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]",
        title=orb.name,
        aspect=DataAspect(),
    )

    function side_ax(row, title_str)
        Axis(fig[row, 2]; title=title_str, xlabel=L"t\ \mathrm{[N.D]}")
    end

    ax_I2 = side_ax(1, L"\hat{I}_2(t)")
    ax_θ2 = side_ax(2, L"\theta_2(t)\ \mathrm{[rad]}")
    ax_I3 = side_ax(3, L"\hat{I}_3(t)")
    ax_θ3 = side_ax(4, L"\theta_3(t)\ \mathrm{[rad]}")

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    lines!(ax_xy,  x,   y;   color=color, linewidth=2)
    lines!(ax_I2,  ts, I2;   color=color, linewidth=1.5)
    lines!(ax_θ2,  ts, θ2;   color=color, linewidth=1.5)
    lines!(ax_I3,  ts, I3;   color=color, linewidth=1.5)
    lines!(ax_θ3,  ts, θ3;   color=color, linewidth=1.5)

    scatter!(ax_xy, [xL1], [0.0]; color=:red, markersize=12)

    if !isnothing(png_path)
        save(png_path, fig; px_per_unit=2)
        println("Saved: $png_path")
    end

    return fig
end


function plot_orbit_actionangle_static(orb::orbit_struct,
                                        transformations;
                                        png_path = nothing,
                                        color    = :royalblue)

    traj = orb.syntraj
    AA   = orb.AAtraj
    ts   = collect(orb.tspan)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)
    z  = getindex.(traj, 3)

    I2 = getindex.(AA, 3)
    θ2 = mod.(getindex.(AA, 4), 2π)
    I3 = getindex.(AA, 5)
    θ3 = mod.(getindex.(AA, 6), 2π)

    xL1 = transformations.recentering[1]

    set_theme!(figure_settings(:large))
    fig = Figure(size=(1600, 900), figure_padding=(80, 20, 20, 20))

    ax3d = Axis3(fig[1:4, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]", zlabel=L"z [N.D]",
        aspect=:data,
        title=orb.name,
        azimuth   = π/3,
        elevation = 0.15π,
    )

    function side_ax(row, title_str)
        Axis(fig[row, 2]; title=title_str, xlabel=L"t\ \mathrm{[N.D]}")
    end

    ax_I2 = side_ax(1, L"\hat{I}_2(t)")
    ax_θ2 = side_ax(2, L"\theta_2(t)\ \mathrm{[rad]}")
    ax_I3 = side_ax(3, L"\hat{I}_3(t)")
    ax_θ3 = side_ax(4, L"\theta_3(t)\ \mathrm{[rad]}")

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    lines!(ax3d,   x,   y,  z;  color=color, linewidth=2)
    lines!(ax_I2,  ts, I2;      color=color, linewidth=1.5)
    lines!(ax_θ2,  ts, θ2;      color=color, linewidth=1.5)
    lines!(ax_I3,  ts, I3;      color=color, linewidth=1.5)
    lines!(ax_θ3,  ts, θ3;      color=color, linewidth=1.5)

    scatter!(ax3d, [xL1], [0.0], [0.0]; color=:red, markersize=12)

    if !isnothing(png_path)
        save(png_path, fig; px_per_unit=2)
        println("Saved: $png_path")
    end

    return fig
end


## ── [1/xx] Load data ──────────────────────────────────────────────────────────
println("Loading data...")
mu = 0.0121505406935702
L1_data_path = "data/Earth-Moon/Resonant/L1/Order16.jld2"
L2_data_path = "data/Earth-Moon/Resonant/L2/Order16.jld2"

transformations = load(L1_data_path, "transformations")
hamiltonians = load(L1_data_path, "hamiltonians")
derivatives = load(L1_data_path, "derivatives")
eoms_NF = derivatives.NF_flow
eoms_AA = derivatives.AA_flow

# transformations = load(L2_data_path, "transformations")
# hamiltonians = load(L2_data_path, "hamiltonians")
# derivatives = load(L2_data_path, "derivatives")
# eoms_NF = derivatives.NF_flow
# eoms_AA = derivatives.AA_flow

println("    L1: $L1_data_path")
println("    L2: $L2_data_path")

## Initialize the halo for a given energy level
println("Initialize halo...")
H = 0.95
θ2dot = eoms_AA[4]
θ3dot = eoms_AA[6]
function residual!(F, x, H_given)
    I2, I3 = x[1], x[2]
    state = SVector(0, 0, I2, pi/2, I3, 0)
    F[1] = evaluatepolynomial(θ2dot, state)
    F[2] = evaluatepolynomial(hamiltonians.action_angle, state)-H_given
end

I_guess = [0.4, 0.4]
prob = NonlinearProblem(residual!, I_guess, H)
sol = solve(prob, NewtonRaphson(); reltol = 1e-12, abstol = 1e-12)

I2_halo, I3_halo = sol.u


lyap_AA = SVector(0, 0, 0.4, 0, 0.4, 0)
lyap_period = 2pi/(evaluatepolynomial(θ2dot,lyap_AA)+evaluatepolynomial(θ3dot,lyap_AA))
halo_AA = SVector(0, 0, I2_halo, pi/2, I3_halo, 0)
halo_period = 2pi/evaluatepolynomial(θ3dot,halo_AA)

quasihalo_N_AA = SVector(0, 0, I2_halo + 0.1, pi/2, I3_halo, 0)
lissajous_AA = SVector(0, 0, I2_halo, 0, I3_halo, 0)

## Propagate orbits
tspans = LinRange(0,100,2000)
tspan_halo = LinRange(0,halo_period,300)
tspan_lyap = LinRange(0,lyap_period,300)

lyap_AAtraj = propagate(eoms_polynomial,lyap_AA,tspan_lyap,eoms_AA)
println("Lyapunov propagated")
halo_AAtraj = propagate(eoms_polynomial,halo_AA,tspan_halo,eoms_AA)
println("Halo propagated")

quasihalo_N_AAtraj = propagate(eoms_polynomial,quasihalo_N_AA,tspans,eoms_AA)
println("QHalo 1 propagated")
lissajous_AAtraj = propagate(eoms_polynomial,lissajous_AA ,tspans,eoms_AA)
println("Lissajous propagated")

## Store the orbits
lyap_orbit = orbit_struct("Lyapunov orbit",tspan_lyap,actionangle2synodic(lyap_AAtraj,transformations),actionangle2normalform(lyap_AAtraj,transformations),lyap_AAtraj)
println("Lyapunov stored")
halo_orbit = orbit_struct("Northern halo orbit",tspan_halo,actionangle2synodic(halo_AAtraj,transformations),actionangle2normalform(halo_AAtraj,transformations),halo_AAtraj)
println("Halo stored")

quasihalo1_orbit = orbit_struct("Northern quasihalo orbit",tspans,actionangle2synodic(quasihalo_N_AAtraj,transformations),actionangle2normalform(quasihalo_N_AAtraj,transformations),quasihalo_N_AAtraj)
println("QHalo 1 stored")
lissajous_orbit = orbit_struct("Lissajous orbit",tspans,actionangle2synodic(lissajous_AAtraj,transformations),actionangle2normalform(lissajous_AAtraj,transformations),lissajous_AAtraj)
println("Lissajous stored")

## Plot and save static PNGs
plot_orbit_phasespace_2d_static(lyap_orbit,transformations; png_path="script/nf_figs/lyapphase_static.png", color=:royalblue)
plot_orbit_actionangle_2d_static(lyap_orbit,transformations; png_path="script/nf_figs/lyapAA_static.png", color=:royalblue)
plot_orbit_phasespace_static(halo_orbit,transformations; png_path="script/nf_figs/halophase_static.png", color=:royalblue)
plot_orbit_actionangle_static(halo_orbit,transformations; png_path="script/nf_figs/haloAA_static.png", color=:royalblue)

plot_orbit_phasespace_static(quasihalo1_orbit,transformations; png_path="script/nf_figs/quasihalo1phase_static.png", color=:royalblue)
plot_orbit_actionangle_static(quasihalo1_orbit,transformations; png_path="script/nf_figs/quasihalo1AA_static.png", color=:royalblue)
plot_orbit_phasespace_static(lissajous_orbit,transformations; png_path="script/nf_figs/lissajousphase_static.png", color=:royalblue)
plot_orbit_actionangle_static(lissajous_orbit,transformations; png_path="script/nf_figs/lissajousAA_static.png", color=:royalblue)


## Synodic halo
halo_syn0 = SVector(0.829990394066360,
0,
0.114000000000000,
0,
0.229330078260162,
0)
synhalo_period = 2.7873
synhalo_syntraj = propagate(eoms_CR3BP,halo_syn0,LinRange(0,3*synhalo_period,1000),mu)
synhalo_NFtraj = synodic2normalform(synhalo_syntraj,transformations)
synhalo_AAtraj = normalform2actionangle(synhalo_NFtraj,transformations)
synhalo_orbit = orbit_struct("Northern halo orbit",LinRange(0,3*synhalo_period,1000),synhalo_syntraj,synhalo_NFtraj,synhalo_AAtraj)
plot_orbit_phasespace_static(synhalo_orbit,transformations; png_path="script/nf_figs/synhalophase_static.png", color=:royalblue)
plot_orbit_actionangle_static(synhalo_orbit,transformations; png_path="script/nf_figs/synhaloAA_static.png", color=:royalblue)
