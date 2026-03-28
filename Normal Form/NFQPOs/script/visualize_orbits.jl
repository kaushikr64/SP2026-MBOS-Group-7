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

## Plot stuff
plot_orbit_phasespace_2d(lyap_orbit,transformations; mp4_path="script/nf_figs/lyapphase.mp4", color=:royalblue,fps = 40)
plot_orbit_actionangle_2d(lyap_orbit,transformations; mp4_path="script/nf_figs/lyapAA.mp4", color=:royalblue,fps = 40)
plot_orbit_phasespace(halo_orbit,transformations; mp4_path="script/nf_figs/halophase.mp4", color=:royalblue, fps = 40)
plot_orbit_actionangle(halo_orbit,transformations; mp4_path="script/nf_figs/haloAA.mp4", color=:royalblue,fps = 40)

plot_orbit_phasespace(quasihalo1_orbit,transformations; mp4_path="script/nf_figs/quasihalo1phase.mp4", color=:royalblue,fps = 90)
plot_orbit_actionangle(quasihalo1_orbit,transformations; mp4_path="script/nf_figs/quasihalo1AA.mp4", color=:royalblue,fps = 90)
plot_orbit_phasespace(lissajous_orbit,transformations; mp4_path="script/nf_figs/lissajousphase.mp4", color=:royalblue,fps = 90)
plot_orbit_actionangle(lissajous_orbit,transformations; mp4_path="script/nf_figs/lissajousAA.mp4", color=:royalblue,fps = 90)


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
plot_orbit_phasespace(synhalo_orbit,transformations; mp4_path="script/nf_figs/synhalophase.mp4", color=:royalblue, fps = 40)
plot_orbit_actionangle(synhalo_orbit,transformations; mp4_path="script/nf_figs/synhaloAA.mp4", color=:royalblue,fps = 40)
