module CR3BPTools

using StaticArrays, DifferentialEquations, LinearAlgebra

include("propagation/eoms.jl")
include("propagation/propagator.jl")

export eoms_CR3BP
export propagate, init_integrator

end # module CR3BPTools

