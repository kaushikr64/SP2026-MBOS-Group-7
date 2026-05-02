using CR3BPTools,
    NFCR3BP,
    StaticArrays,
    LinearAlgebra,
    DifferentialEquations,
    JLD2,
    GLMakie,
    LaTeXStrings

include("../util/plotting/plotting.jl")

## ── [1] Load data ─────────────────────────────────────────────────────────────
println("[1] Loading data...")

L1_data_path = "./data/Earth-Moon/RCM/L1/Order16.jld2"
transformations = load(L1_data_path, "transformations")
derivatives     = load(L1_data_path, "derivatives")
eoms_NF         = derivatives.NF_flow

@load "script/RCM-poincare.jld2" poincare_points

## ── [2] Find non-empty trajectories ──────────────────────────────────────────
println("[2] Enumerating trajectories...")

nonempty_traj_idxs = [i for (i, pts) in enumerate(poincare_points) if !isempty(pts)]

println("    Non-empty trajectories: $(length(nonempty_traj_idxs))")

## ── [3] User selection ────────────────────────────────────────────────────────
selected_idx    = 31        # ← change this (1 … $(length(nonempty_traj_idxs)))
highlight_color = RGBf(0.8, 0.6, 0.1)

sel_traj_idx    = nonempty_traj_idxs[selected_idx]
sel_NF          = poincare_points[sel_traj_idx][1]

## ── [4] Integrate orbit from selected crossing ────────────────────────────────
println("[4] Integrating orbit from crossing $selected_idx...")

orbit_NF  = propagate(eoms_polynomial, sel_NF, LinRange(0, 200, 10000), eoms_NF)
orbit_syn = normalform2synodic(orbit_NF, transformations)
filter!(s -> norm(s[1:3]) <= 2, orbit_syn)

## ── [5] Plot ──────────────────────────────────────────────────────────────────
println("[5] Plotting...")

xL1 = transformations.recentering[1]

set_theme!(figure_settings(:large))

# ── Figure 1: 2D Poincaré map ─────────────────────────────────────────────────
fig_map = Figure()
ax_map = Axis(
    fig_map[1, 1];
    xlabel = L"$x$ [n.d]",
    ylabel = L"$y$ [n.d]",
    aspect = DataAspect(),
    title  = L"Poincaré map — trajectory $%$(sel_traj_idx)$",
)

for (i, pts) in enumerate(poincare_points)
    isempty(pts) && continue
    pts_syn = normalform2synodic(pts, transformations)
    pts_syn = filter(s -> norm(s[1:3]) <= 2, pts_syn)
    isempty(pts_syn) && continue
    color = i == sel_traj_idx ? highlight_color : :black
    scatter!(ax_map, getindex.(pts_syn, 1), getindex.(pts_syn, 2);
        color = color, markersize = 5)
end

display(fig_map)

# ── Figure 2: 3D orbit + full Poincaré map ────────────────────────────────────
fig_3d = Figure(figure_padding = (100, 80, 60, 40))
ax_3d = Axis3(
    fig_3d[1, 1];
    xlabel             = L"$x$ [n.d]",
    ylabel             = L"$y$ [n.d]",
    zlabel             = L"$z$ [n.d]",
    aspect             = :data,
    title              = L"Orbit from trajectory $%$(sel_traj_idx)$",
    azimuth            = π / 4,
    elevation          = 0.15π,
    xautolimitmargin   = (0.1, 0.1),
    yautolimitmargin   = (0.1, 0.1),
    zautolimitmargin   = (0.15, 0.15),
)

for (i, pts) in enumerate(poincare_points)
    isempty(pts) && continue
    pts_syn = normalform2synodic(pts, transformations)
    pts_syn = filter(s -> norm(s[1:3]) <= 2, pts_syn)
    isempty(pts_syn) && continue
    color = i == sel_traj_idx ? highlight_color : (:black, 0.4)
    scatter!(ax_3d,
        getindex.(pts_syn, 1), getindex.(pts_syn, 2), getindex.(pts_syn, 3);
        color = color, markersize = i == sel_traj_idx ? 8 : 4, transparency = true)
end

lines!(ax_3d,
    getindex.(orbit_syn, 1), getindex.(orbit_syn, 2), getindex.(orbit_syn, 3);
    color = (highlight_color, 0.8), linewidth = 1.5)

display(fig_3d)


# [RGBf(0.25, 0.45, 0.75), RGBf(0.75, 0.28, 0.1), RGBf(0.8, 0.6, 0.1)]