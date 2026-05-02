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
    NonlinearSolve,
    ProgressMeter

include("../util/plotting/plotting.jl")

## ── [1/xx] Load data ──────────────────────────────────────────────────────────
println("[1] Loading data...")

L1_data_path = "data/Earth-Moon/Resonant/L1/Order16.jld2"
mu = 0.012150585609624

transformations = load(L1_data_path, "transformations")
hamiltonians = load(L1_data_path, "hamiltonians")
derivatives = load(L1_data_path, "derivatives")
eoms_NF = derivatives.NF_flow
eoms_AA = derivatives.AA_flow

println("    L1: $L1_data_path")


## Initialize the halo for a given energy level
H = LinRange(0.4,10,1000)
θ2dot = eoms_AA[4]
function residual!(F, x, H_given)
    I2, I3 = x[1], x[2]
    state = SVector(0, 0, I2, pi/2, I3, 0)
    F[1] = evaluatepolynomial(θ2dot, state)
    F[2] = evaluatepolynomial(hamiltonians.action_angle, state)-H_given
end

halo_AA_list = []
JCs = []
JC_calc = x -> (x[1]^2 + x[2]^2 + 2*(1-mu)/norm([x[1]+mu,x[2],x[3]]) + 2*mu/norm([x[1]+mu,x[2],x[3]])-dot([x[4],x[5],x[6]],[x[4],x[5],x[6]]))

I_guess = [0.4, 0.4]
for i in eachindex(H)
    prob = NonlinearProblem(residual!, I_guess, H[i])
    sol = solve(prob, NewtonRaphson(); reltol = 1e-12, abstol = 1e-12)
    I2_halo, I3_halo = sol.u
    halo_AA = SVector(0, 0, I2_halo, pi/2, I3_halo, 0)
    push!(halo_AA_list,halo_AA)
    halo_syn = actionangle2synodic(halo_AA,transformations)
    JC = JC_calc(halo_syn)
    push!(JCs,JC)
end


## ── [2/xx] Propagate halos ────────────────────────────────────────────────
println("[2] Propagating halos...")
θ3dot = eoms_AA[6]
stride = 10  # plot every 10th orbit (~100 total)
syn_trajs = []
jc_subset = Float64[]

@showprogress for i in 1:stride:length(halo_AA_list)
    halo_AA = halo_AA_list[i]
    period = 2π / evaluatepolynomial(θ3dot, halo_AA)
    tspan = LinRange(0, period, 200)
    AAtraj = propagate(eoms_polynomial, halo_AA, tspan, eoms_AA)
    syntraj = actionangle2synodic(AAtraj, transformations)
    push!(syn_trajs, syntraj)
    push!(jc_subset, JCs[i])
end

## ── [3/xx] Plot halo family ───────────────────────────────────────────────
println("[3] Plotting halo family...")
fig = Figure(size=(900, 700))
ax = Axis3(fig[1, 1];
    xlabel="x [N.D]", ylabel="y [N.D]", zlabel="z [N.D]",
    title="Halo Family — Earth-Moon L1",
    aspect=:data,
)

jc_min, jc_max = extrema(jc_subset)
cmap = :viridis
for (traj, jc) in zip(syn_trajs, jc_subset)
    c = (jc - jc_min) / (jc_max - jc_min)
    x = [s[1] for s in traj]
    y = [s[2] for s in traj]
    z = [s[3] for s in traj]
    lines!(ax, x, y, z; color=c, colorrange=(0, 1), colormap=cmap, linewidth=0.8)
end
Colorbar(fig[1, 2]; colormap=cmap, limits=(jc_min, jc_max), label="Jacobi Constant")

# save("script/nf_figs/halo_family.png", fig)
# println("Saved: script/nf_figs/halo_family.png")

