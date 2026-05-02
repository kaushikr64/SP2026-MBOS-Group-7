using CR3BPTools,
    NFCR3BP,
    StaticArrays,
    LinearAlgebra,
    DifferentialEquations,
    Roots,
    JLD2,
    MathOptInterface,
    GLMakie,
    LaTeXStrings,
    MATLAB,
    NonlinearSolve,
    ProgressMeter,
    Printf

include("../util/plotting/plotting.jl")

## ── [1/xx] Load data ──────────────────────────────────────────────────────────
println("[1] Loading data...")

L1_data_path = "./data/Earth-Moon/RCM/L1/Order16.jld2"

transformations = load(L1_data_path, "transformations")
hamiltonians = load(L1_data_path, "hamiltonians")
derivatives = load(L1_data_path, "derivatives")
eoms_NF = derivatives.NF_flow

mat"load('script/L1QuasiHaloFam.mat')"
mat"addpath('script')"

fv_hst = mat"L1QuasiHaloFam.fv_hst"   # 20×3×367 Float64
JC_hst = mat"L1QuasiHaloFam.JC_hst"   # 20×1 Float64
X0_hst = mat"L1QuasiHaloFam.X0_hst"   # 20×6 Float64

## ── [2/xx] Define structures ──────────────────────────────────────────────────

# --- Frequency data (from MATLAB) ---

struct GivenFreqs
    fv::Vector{Float64}          # 367 frequencies
end

struct GivenEnergyFreqs
    JC::Float64
    x0::SVector{6,Float64}
    orbits::Vector{GivenFreqs}   # 3 orbits
end

struct AllFreqs
    energies::Vector{GivenEnergyFreqs}  # 20 energy levels
end

function AllFreqs(fv_hst, JC_hst, X0_hst)
    n_energy = size(fv_hst, 1)  # 20
    n_orbs   = size(fv_hst, 2)  # 3
    energies = [
        GivenEnergyFreqs(
            JC_hst[i],
            SVector{6,Float64}(X0_hst[i, :]),
            [GivenFreqs(fv_hst[i, j, :]) for j in 1:n_orbs]
        )
        for i in 1:n_energy
    ]
    return AllFreqs(energies)
end

# --- Computed map sections ---

struct OrbitSections
    poincare_bg::Vector{SVector{6,Float64}}   # Poincaré map crossings
    poincare_orbits::Vector{SVector{6,Float64}}   # Poincaré map crossings
    kkg::Vector{SVector{6,Float64}}        # KKG section points
end

struct EnergyOrbits
    JC::Float64
    x0::SVector{6,Float64}
    orbits::Vector{OrbitSections}          # 3 orbits
end

struct AllOrbits
    energies::Vector{EnergyOrbits}         # 20 energy levels
end

function AllOrbits(all_freqs::AllFreqs, energy_indices)
    energies = [
        EnergyOrbits(
            all_freqs.energies[i].JC,
            all_freqs.energies[i].x0,
            [OrbitSections(SVector{6,Float64}[], SVector{6,Float64}[], SVector{6,Float64}[]) for _ in all_freqs.energies[i].orbits]
        )
        for i in energy_indices
    ]
    return AllOrbits(energies)
end

function compute_poincare_map(H, nf0, hamiltonians, eoms_NF; numcrosses=200, n_pts=5, q2_range=0.4)
    q1, p1, q3 = 0.0, 0.0, 0.0
    q2_axis_pts = LinRange(nf0[3], nf0[3] - q2_range, n_pts)

    p3_calc_fn = (q2, p3) ->
        evaluatepolynomial(hamiltonians.real_normal, SVector(q1, p1, q2, 0.0, q3, p3)) - H

    poincare_points = [SVector{6,Float64}[] for _ in eachindex(q2_axis_pts)]

    prog = Progress(
        length(q2_axis_pts) * numcrosses;
        desc = "Computing Poincaré map...",
        output = stderr,
    )

    condition(u, t, integrator) = u[5]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(
        condition, affect!, nothing;
        rootfind = SciMLBase.RightRootFind,
        save_positions = (false, false),
    )

    IC0   = SVector(q1, p1, q2_axis_pts[1], 0.0, q3, 0.0)
    tspan = (0.0, 10.0)
    prob0 = ODEProblem(eoms_polynomial, IC0, tspan, eoms_NF)

    Threads.@threads for i in eachindex(q2_axis_pts)
        q2 = q2_axis_pts[i]
        p3 = find_zero(p3 -> p3_calc_fn(q2, p3), 0.1)
        IC = SVector(q1, p1, q2, 0.0, q3, p3)

        crossings = sizehint!(SVector{6,Float64}[], numcrosses)
        for _ in 1:numcrosses
            prob = remake(prob0; u0 = IC, tspan = tspan)
            sol  = solve(prob, Vern9();
                callback = cb, save_everystep = false, dense = false,
                abstol = 1e-9, reltol = 1e-9)
            if norm(sol.u[end]) < 10
                IC = sol.u[end]
                push!(crossings, IC)
                next!(prog)
            else
                break
            end
        end
        poincare_points[i] = crossings
    end
    finish!(prog)
    return poincare_points
end

function compute_poincare_ring(nf0, eoms_NF; numcrosses=200)
    condition(u, t, integrator) = u[5]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(
        condition, affect!, nothing;
        rootfind = SciMLBase.RightRootFind,
        save_positions = (false, false),
    )

    tspan = (0.0, 10.0)
    prob0 = ODEProblem(eoms_polynomial, nf0, tspan, eoms_NF)

    IC = SVector(0,0,nf0[3:6]...)
    crossings = sizehint!(SVector{6,Float64}[], numcrosses)
    for _ in 1:numcrosses
        prob = remake(prob0; u0 = IC, tspan = tspan)
        sol  = solve(prob, Vern9();
            callback = cb, save_everystep = false, dense = false,
            abstol = 1e-9, reltol = 1e-9)
        if norm(sol.u[end]) < 10
            IC = sol.u[end]
            push!(crossings, IC)
        else
            break
        end
    end
    return crossings
end

## ── [3/xx] Build structures ───────────────────────────────────────────────────
println("[3] Building structures...")

n_energy_levels = 10   # set to e.g. 5 to use only first 5; nothing = all 20

all_freqs      = AllFreqs(fv_hst, JC_hst, X0_hst)
energy_indices = isnothing(n_energy_levels) ? eachindex(all_freqs.energies) : 1:n_energy_levels
all_orbits     = AllOrbits(all_freqs, energy_indices)

## ── [4/xx] Compute maps ───────────────────────────────────────────────────────
println("[4] Computing maps...")

@showprogress for (local_i, orig_i) in enumerate(energy_indices)
    nf0 = synodic2normalform(all_orbits.energies[local_i].x0, transformations)
    H = evaluatepolynomial(hamiltonians.real_normal,nf0)
    # poincare_points = compute_poincare_map(H, nf0, hamiltonians, eoms_NF)

    for j in 1:3
        println(j)
        # for pts in poincare_points
        #     append!(all_orbits.energies[local_i].orbits[j].poincare_bg, normalform2synodic(pts,transformations))
        # end
        kkg_mat = mat"KKG_fv_2_SecHst(L1QuasiHaloFam.fv_hst($orig_i, $j, :))"
        for k in axes(kkg_mat, 2)
            push!(all_orbits.energies[local_i].orbits[j].kkg, SVector{6,Float64}(kkg_mat[:, k]))
        end
        println(j)
        x0_nf = synodic2normalform(SVector{6,Float64}(kkg_mat[:, 1]), transformations)
        ring   = compute_poincare_ring(x0_nf, eoms_NF)
        append!(all_orbits.energies[local_i].orbits[j].poincare_orbits,
                normalform2synodic(ring, transformations))
    end
end

## ── [5/xx] Plot KKG sections ──────────────────────────────────────────────────
println("[5] Plotting KKG sections...")

# Pre-extract all synodic data per energy level
kkg_syn = [
    [filter(s -> norm(s[1:3]) <= 2, all_orbits.energies[i].orbits[j].kkg)
     for j in 1:3]
    for i in eachindex(all_orbits.energies)
]
bg_syn = [
    filter(s -> norm(s[1:3]) <= 2, all_orbits.energies[i].orbits[1].poincare_bg)
    for i in eachindex(all_orbits.energies)
]
orb_syn = [
    [filter(s -> norm(s[1:3]) <= 2, all_orbits.energies[i].orbits[j].poincare_orbits)
     for j in 1:3]
    for i in eachindex(all_orbits.energies)
]

# Global axis limits across all energy levels
all_x = [s[1] for e in orb_syn for traj in e for s in traj]
all_y = [s[2] for e in orb_syn for traj in e for s in traj]
pad(lo, hi) = (p = 0.07 * (hi - lo); (lo - p, hi + p))
xrng = pad(extrema(all_x)...)
yrng = pad(extrema(all_y)...)

orbit_colors = [:royalblue, :crimson, :forestgreen]

legend_entries = [
    [MarkerElement(color=c, marker='o', markersize=10) for c in orbit_colors]...,
    MarkerElement(color=:black, marker=:circle, markersize=14),
]
legend_labels = ["Quasihalo 1", "Quasihalo 2", "Quasihalo 3", "Halo orbit"]

function make_fig_ax()
    fig = Figure(size=(900, 750))
    ax  = Axis(fig[1, 1];
        xlabel = L"$x$ [n.d]",
        ylabel = L"$y$ [n.d]",
        aspect = DataAspect(),
    )
    Legend(fig[1, 2], legend_entries, legend_labels)
    return fig, ax
end

function draw_frame!(ax, i)
    empty!(ax)
    xlims!(ax, xrng...)
    ylims!(ax, yrng...)
    JC_str = @sprintf("%.4f", all_orbits.energies[i].JC)
    ax.title[] = L"$JC = %$(JC_str)$"

    # Background Poincaré map (disabled)
    # bg = bg_syn[i]
    # if !isempty(bg)
    #     scatter!(ax, getindex.(bg, 1), getindex.(bg, 2);
    #              color=(:gray60, 0.4), markersize=3)
    # end

    for j in 1:3
        orb = orb_syn[i][j]
        !isempty(orb) && scatter!(ax, getindex.(orb, 1), getindex.(orb, 2);
                                   color=orbit_colors[j], markersize=10, marker='o')
        kkg = kkg_syn[i][j]
        !isempty(kkg) && lines!(ax, getindex.(kkg, 1), getindex.(kkg, 2);
                                  color=orbit_colors[j], linewidth=1.5)
    end
    x0 = all_orbits.energies[i].x0
    scatter!(ax, [x0[1]], [x0[2]]; color=:black, marker=:circle, markersize=14)
end

set_theme!(figure_settings(:large))

# Open each energy level as a separate interactive figure
for i in eachindex(all_orbits.energies)
    fig_i, ax_i = make_fig_ax()
    draw_frame!(ax_i, i)
    display(GLMakie.Screen(), fig_i)
end

# Save animation
fig, ax = make_fig_ax()
record(fig, "script/step3_figs/kkg_sections.mp4", eachindex(all_orbits.energies); framerate=1) do i
    draw_frame!(ax, i)
end
println("Saved: script/step3_figs/kkg_sections.mp4")
