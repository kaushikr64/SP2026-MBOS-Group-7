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
    ProgressMeter

include("../util/plotting/plotting.jl")

## ── [1/xx] Load data ──────────────────────────────────────────────────────────
println("[1] Loading data...")

L1_data_path = "./data/Earth-Moon/RCM/L1/Order16.jld2"

transformations = load(L1_data_path, "transformations")
hamiltonians = load(L1_data_path, "hamiltonians")
derivatives = load(L1_data_path, "derivatives")
eoms_NF = derivatives.NF_flow

# transformations = load(L2_data_path, "transformations")
# hamiltonians = load(L2_data_path, "hamiltonians")
# derivatives = load(L2_data_path, "derivatives")
# eoms_NF = derivatives.NF_flow
# eoms_AA = derivatives.AA_flow

println("    L1: $L1_data_path")

## Initialize the halo for a given energy level
println("[2] Initialize halo...")
H = 0.9


## Start setting up the poincaré map
println("[3] Set up Poincaré map...")
numcrosses = 300;
q1 = 0
p1 = 0
q2_axis_pts = LinRange(-0.85, 0.85, 51)
p2_axis_pts = LinRange(-0.25, 0, 3)
q3 = 0
p3_calc_fn_q2 =
    (q2, p3) ->
        evaluatepolynomial(hamiltonians.real_normal, SVector(q1, p1, q2, 0, q3, p3))-H
p3_calc_fn_p2 =
    (p2, p3) ->
        evaluatepolynomial(hamiltonians.real_normal, SVector(q1, p1, 0, p2, q3, p3))-H        

# Initialize the poincare map points
poincare_points = [SVector{6,Float64}[] for _ = 1:length(q2_axis_pts)+length(p2_axis_pts)]


## Run poincare map
println("[4] Running Poincaré map...")
prog = Progress(
    (length(q2_axis_pts)+length(p2_axis_pts)) * numcrosses;
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
IC0 = SVector(q1, p1, q2_axis_pts[1], 0, q3, 0.0)
tspan = (0.0, 10.0)
prob0 = ODEProblem(eoms_polynomial, IC0, tspan, eoms_NF)

Threads.@threads for i in eachindex(q2_axis_pts)
    q2 = q2_axis_pts[i]
    p3 = find_zero(p3 -> p3_calc_fn_q2(q2, p3), 0.1)
    IC = SVector(q1, p1, q2, 0, q3, p3)

    crossings = sizehint!(SVector{6,Float64}[], numcrosses)

    for j = 1:numcrosses
        prob = remake(prob0; u0 = IC, tspan = tspan)
        sol = solve(
            prob,
            Vern9();
            callback = cb,
            save_everystep = false,
            dense = false,
            abstol = 1e-9,
            reltol = 1e-9,
        )
        if norm(sol.u[end])<10
            IC = sol.u[end]
            push!(crossings, IC)
            next!(prog)
        else
            break
        end
    end
    poincare_points[i] = crossings
end

Threads.@threads for i in (eachindex(p2_axis_pts))
    p2 = p2_axis_pts[i]
    p3 = find_zero(p3 -> p3_calc_fn_p2(p2, p3), 0.1)
    IC = SVector(q1, p1, 0, p2, q3, p3)

    crossings = sizehint!(SVector{6,Float64}[], numcrosses)

    for j = 1:numcrosses
        prob = remake(prob0; u0 = IC, tspan = tspan)
        sol = solve(
            prob,
            Vern9();
            callback = cb,
            save_everystep = false,
            dense = false,
            abstol = 1e-9,
            reltol = 1e-9,
        )
        if norm(sol.u[end])<10
            IC = sol.u[end]
            push!(crossings, IC)
            next!(prog)
        else
            break
        end
    end
    poincare_points[i + length(q2_axis_pts)] = crossings
end
finish!(prog)

## Plot the poincaré map
@load "script/RCM-poincare.jld2" poincare_points
JC_calc = x -> (x[1]^2 + x[2]^2 + 2*(1-mu)/norm([x[1]+mu,x[2],x[3]]) + 2*mu/norm([x[1]+mu,x[2],x[3]])-dot([x[4],x[5],x[6]],[x[4],x[5],x[6]]))
JC = JC_calc(normalform2synodic(poincare_points[findfirst(!isempty, poincare_points)][1],transformations))

set_theme!(figure_settings(:large))
fig_poincare_nf = Figure()
ax = Axis(
    fig_poincare_nf[1, 1];
    xlabel = L"$x$ [n.d]",
    ylabel = L"$y$ [n.d]",
    aspect = DataAspect(),
    title = L"$JC = %$(round(JC,digits = 3))$"
)

for pts in poincare_points
    isempty(pts) && continue
    pts_syn = normalform2synodic(pts, transformations)
    pts_syn = filter(s -> norm(s[1:3]) <= 2, pts_syn)
    isempty(pts_syn) && continue
    scatter!(ax, getindex.(pts_syn, 1), getindex.(pts_syn, 2); color=:black, markersize=5)
end


display(fig_poincare_nf)

@save "script/RCM-poincare.jld2" poincare_points


