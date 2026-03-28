using CR3BPTools,
    NFCR3BP,
    StaticArrays,
    LinearAlgebra,
    DifferentialEquations,
    Roots,
    JLD2,
    MathOptInterface,
    Snopt,
    SNOW,
    GLMakie,
    ForwardDiff,
    LaTeXStrings,
    MATLAB,
    NonlinearSolve,
    ProgressMeter

include("../util/plotting/plotting.jl")

## ── [1/xx] Load data ──────────────────────────────────────────────────────────
println("[1] Loading data...")

L1_data_path = "../data/Earth-Moon/Resonant/L1/Order16.jld2"
L2_data_path = "../data/Earth-Moon/Resonant/L2/Order16.jld2"

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
println("[2] Initialize halo...")
H = 0.85
θ2dot = eoms_AA[4]
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

halo_AA = SVector(0, 0, I2_halo, pi/2, I3_halo, 0)
halo_NF = actionangle2normalform(halo_AA, transformations)

## Start setting up the poincaré map
println("[3] Set up Poincaré map...")
numcrosses = 2000;
q1 = 0
p1 = 0
q2_axis_pts = LinRange(-halo_NF[4], halo_NF[4], 51);
p2 = 0
q3 = 0
p3_calc_fn =
    (q2, p3) ->
        evaluatepolynomial(hamiltonians.real_normal, SVector(q1, p1, q2, p2, q3, p3))-H

# Initialize the poincare map points
poincare_points = [SVector{6,Float64}[] for _ = 1:length(q2_axis_pts)]




## Run poincare map
println("[4] Running Poincaré map...")
prog = Progress(
    length(q2_axis_pts) * numcrosses;
    desc = "Computing Poincaré map...",
    output = stderr,
)

# Define callback once, outside all loops
condition(u, t, integrator) = u[5]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(
    condition,
    affect!,
    nothing;
    rootfind = SciMLBase.RightRootFind,
    save_positions = (false, false),
)

# Build a template problem once
IC0 = SVector(q1, p1, q2_axis_pts[1], p2, q3, 0.0)
tspan = (0.0, 10.0)
prob0 = ODEProblem(eoms_polynomial, IC0, tspan, eoms_NF)

Threads.@threads for i in eachindex(q2_axis_pts)
    q2 = q2_axis_pts[i]
    p3 = find_zero(p3 -> p3_calc_fn(q2, p3), 0.1)
    IC = SVector(q1, p1, q2, p2, q3, p3)
    
    crossings = sizehint!(SVector{6,Float64}[], numcrosses)

    for j = 1:numcrosses
        prob = remake(prob0; u0 = IC, tspan = tspan)
        sol = solve(
            prob,
            Vern9();
            callback = cb,
            save_everystep = false,
            dense = false,
            abstol = 1e-6,   
            reltol = 1e-6,
        )
        IC = sol.u[end]
        push!(crossings, IC)
        next!(prog)
    end
    poincare_points[i] = crossings
end
finish!(prog)

## Plot the poincaré map
set_theme!(figure_settings(:large))
fig_poincare_nf = Figure()
ax = Axis(
    fig_poincare_nf[1, 1];
    xlabel = L"$\bar{q}_2$ [n.d]",
    ylabel = L"$\bar{p}_2$ [n.d]",
    # xticks = (range(0, 2π, length=5),
    #          ["0", "π/2", "π", "3π/2", "2π"])
)

cmap = cgrad(:viridis)
I3_vals = [pts[1][5] for pts in poincare_points if !isempty(pts)]
clims = (minimum(I3_vals), maximum(I3_vals))


for pts in poincare_points
    I3 = pts[1][5]
    c_ratio = (I3-clims[1])/(clims[2]-clims[1])
    scatter!(ax, getindex.(pts, 3), getindex.(pts, 4), color = cmap[c_ratio])
end

Colorbar(fig_poincare_nf[1, 2]; colormap = :viridis, limits = clims, label = L"\hat{I}_3")
display(fig_poincare_nf)


## Plot the poincaré map
set_theme!(figure_settings(:large))
@load "script/Analysis/data/poincare.jld2" poincare_aapoints poincare_points

poincare_aapoints = [
    normalform2actionangle(poincare_point, transformations) for
    poincare_point in poincare_points
]

fig_poincare_aa = Figure(size = (850,550),figure_padding=(10, 60, 10, 10))
ax = Axis(
    fig_poincare_aa[1, 1];
    xlabel = L"$\theta_2$ [rad]",
    ylabel = L"$\hat{I}_2$ [n.d]",
    xticks = (range(0, 2π, length = 5), ["0", "π/2", "π", "3π/2", "2π"]),
    title = L"H = 0.95",
)

cmap = cgrad(:viridis)
I3_vals = [pts[1][5] for pts in poincare_aapoints if !isempty(pts)]
clims = (minimum(I3_vals), maximum(I3_vals))

# Plot everything but vertical
for i in eachindex(poincare_aapoints)
    pts = poincare_aapoints[i]
    if i != 26
        if i == 1
            scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :blue, markersize = 12,label = "Northern halo")
        elseif i in 2:9
            if i == 2
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :steelblue, markersize = 7,label = "Northern quaishalo")
            else
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :steelblue, markersize = 7)
            end
        elseif i in 43:50
            if i == 43
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :orange, markersize = 7,label = "Southern quaishalo")
            else
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :orange, markersize = 7)
            end
        elseif i == 51
            scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color = :red, markersize = 12,label = "Southern halo")
        else
            if i == 10
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color=(:gray60), markersize = 7,label = "Lissajous")
            else
                scatter!(ax, getindex.(pts, 4), getindex.(pts, 3), color=(:gray60), markersize = 7)
            end
        end
    end
end
# Plot verticals at the end
let 
    pts = poincare_aapoints[26]
    scatter!(ax, LinRange(0,2π,length(pts)), getindex.(pts, 3), color = :purple, markersize = 14,label = "Vertical")
end

Legend(fig_poincare_aa[2, 1], ax; orientation=:horizontal, nbanks=2)

tightlimits!(ax)
ylims!(ax, nothing, 0.42) 
display(fig_poincare_aa)
save("script/Analysis/fig/poincare_aa.png", fig_poincare_aa; px_per_unit = 5)

@save "script/Analysis/data/poincare.jld2" poincare_aapoints poincare_points

## Plot the poincaré map
set_theme!(figure_settings(:large))
@load "script/Analysis/data/poincare.jld2" poincare_aapoints poincare_points

poincare_synpoints = [
    normalform2synodic(poincare_point, transformations) for
    poincare_point in poincare_points
]

fig_poincare_syn = Figure(size = (850,1050),figure_padding=(10, 60, 10, 10))
ax = Axis(
    fig_poincare_syn[1, 1];
    xlabel = L"$x$ [n.d]",
    ylabel = L"$y$ [n.d]",
    title = L"H = 0.95",
    aspect = DataAspect()
)

cmap = cgrad(:viridis)
I3_vals = [pts[1][5] for pts in poincare_synpoints if !isempty(pts)]
clims = (minimum(I3_vals), maximum(I3_vals))

# Plot everything but vertical
for i in eachindex(poincare_synpoints)
    pts = poincare_synpoints[i]
    if i == 1
        scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :blue, markersize = 12,label = "Northern halo")
    elseif i in 2:9
        if i == 2
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :steelblue, markersize = 7,label = "Northern quaishalo")
        else
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :steelblue, markersize = 7)
        end
    elseif i in 43:50
        if i == 43
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :orange, markersize = 7,label = "Southern quaishalo")
        else
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :orange, markersize = 7)
        end
    elseif i == 51
        scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :red, markersize = 12,label = "Southern halo")
    elseif i == 26
        scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color = :purple, markersize = 12,label = "Vertical")
    else
        if i == 10
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color=(:gray60), markersize = 7,label = "Lissajous")
        else
            scatter!(ax, getindex.(pts, 1), getindex.(pts, 2), color=(:gray60), markersize = 7)
        end
    end
end

# Legend(fig_poincare_syn[1, 2], ax; orientation=:vertical)

display(fig_poincare_syn)
save("script/Analysis/fig/poincare_syn.png", fig_poincare_syn; px_per_unit = 5)

@save "script/Analysis/data/poincare.jld2" poincare_aapoints poincare_synpoints poincare_points
