using CR3BPTools,
    NFCR3BP,
    StaticArrays,
    LinearAlgebra,
    DifferentialEquations,
    Roots,
    JLD2,
    NonlinearSolve,
    GLMakie,
    LaTeXStrings,
    Statistics

include("../util/plotting/plotting.jl")

# ── [1/xx] Load data ──────────────────────────────────────────────────────────
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

## Initialize orbits in the normal form
println("Initialize orbits...")
H = 0.95
θ2dot = eoms_AA[4]
θ3dot = eoms_AA[6]
function halo_residual!(F, x, H_given)
    I2, I3 = x[1], x[2]
    state = SVector(0, 0, I2, pi/2, I3, 0)
    F[1] = evaluatepolynomial(θ2dot, state)
    F[2] = evaluatepolynomial(hamiltonians.action_angle, state)-H_given
end

I_guess = [0.2, 0.2]
prob = NonlinearProblem(halo_residual!, I_guess, H)
sol = solve(prob, NewtonRaphson(); reltol = 1e-12, abstol = 1e-12)

I2_halo, I3_halo = sol.u


halo_AA = SVector(0, 0, I2_halo, pi/2, I3_halo, 0)
halo_period = 2pi/evaluatepolynomial(θ3dot,halo_AA)

quasihalo_N_AA = SVector(0, 0, I2_halo + 0.02, pi/2, I3_halo, 0)
quasihalo_N_Syn = actionangle2synodic(quasihalo_N_AA,transformations)

## Converge?
revs = 10
period_rough = 2*pi/evaluatepolynomial(θ3dot, quasihalo_N_AA)

function qpo_residual!(F, X, I_target)
    T = X[end]
    I3_total = 0
    T_total = 0
    for i in 1:revs
        idx      = 5*(i-1)+1
        idx_next = 5*(i % revs)+1
        
        xi      = X[idx:idx+4]
        xi_next = X[idx_next:idx_next+4]
        
        x0_i    = SVector(xi[1],      0.0, xi[2],      xi[3],      xi[4],      xi[5])
        x0_next = SVector(xi_next[1], 0.0, xi_next[2], xi_next[3], xi_next[4], xi_next[5])
        
        crossing = ContinuousCallback(
            (u, t, integrator) -> u[2],
            terminate!,
            affect_neg! = nothing
        )
        traj    = propagate(eoms_CR3BP, x0_i, (0.0, 2*T), mu; callback=crossing)
        xf      = traj.u[end]
        traj_pts = length(traj)

        # state continuity
        eq = 5*(i-1)
        F[eq+1] = xf[1] - x0_next[1]
        F[eq+2] = xf[3] - x0_next[3]
        F[eq+3] = xf[4] - x0_next[4]
        F[eq+4] = xf[5] - x0_next[5]
        F[eq+5] = xf[6] - x0_next[6]

        # accumulate action integral for this arc
        traj_aa  = synodic2actionangle(traj.u, transformations)
        tdiffs   = [traj.t[k+1] - traj.t[k] for k in 1:traj_pts-1]
        arc_time = traj.t[end] - traj.t[1]
        I3_arc   = sum(getindex.(traj_aa, 5)[1:end-1] .* tdiffs)
        
        
        I3_total += I3_arc
        T_total  += arc_time

    end
    
    F[5*revs+1] = (I3_total / T_total) - I_target
end

# Initial guess: distribute revs points around the QPO
# Use normal form to get synodic state, all start near same point
X_guess = zeros(5*revs + 1)
x0_syn = actionangle2synodic(quasihalo_N_AA, transformations)
for i in 1:revs
    idx = 5*(i-1)+1
    X_guess[idx:idx+4] = [x0_syn[1], x0_syn[3], x0_syn[4], x0_syn[5], x0_syn[6]]
end
X_guess[end] = period_rough

prob_diffcorr = NonlinearProblem(qpo_residual!, X_guess, I3_halo + 0.02)
solqpo = solve(prob_diffcorr, NewtonRaphson(), show_trace=Val(true))

## Verify???
T = solqpo.u[end]
x0 = SVector(solqpo.u[1], 0.0, solqpo.u[2], solqpo.u[3], solqpo.u[4], solqpo.u[5])
tspan = LinRange(0, revs*T, 1000)
syntraj = propagate(eoms_CR3BP, x0, tspan, mu)

set_theme!(figure_settings(:small))
fig = Figure(size=(1600, 800), figure_padding=(80, 20, 20, 20))
ax = Axis3(fig[1, 1];
    xlabel=L"x [N.D]", ylabel=L"y [N.D]", zlabel=L"z [N.D]",
    aspect=:data,
)
lines!(ax, getindex.(syntraj, 1), getindex.(syntraj, 2), getindex.(syntraj, 3))
display(fig)
