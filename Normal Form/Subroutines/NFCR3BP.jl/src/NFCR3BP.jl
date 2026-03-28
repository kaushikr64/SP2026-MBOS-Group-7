module NFCR3BP

using DifferentialEquations, StaticArrays, SparseArrays, LinearAlgebra, ProgressMeter, JLD2

include("polynomials/polynomial_definitions.jl")
include("polynomials/polynomial_indices.jl")
include("polynomials/polynomial_operations.jl")
include("polynomials/polynomial_utilities.jl")

include("physics/eoms.jl")
include("physics/find_gamma.jl")
include("physics/get_hamiltonian_derivatives.jl")

include("transformations/recentering.jl")
include("transformations/diagonalize.jl")
include("transformations/realcomplexify.jl")
include("transformations/lieseries.jl")
include("transformations/actionangles.jl")
include("transformations/coordinate_changes.jl")

include("generate_nf/generate_NF.jl")
include("generate_nf/compute_store_NF.jl")

export AbstractHomogenousPolynomial,
    AbstractActionAnglePolynomial,
    RealHomogenousPolynomial,
    ComplexHomogenousPolynomial,
    MixedDegreePolynomial,
    ResonantActionAnglePolynomial,
    ResonantActionAnglePolynomialDegree,
    CombinePolynomials

export PSI,
    MIDX3,
    MIDX6,
    MAXDEGREE,
    get_multiindex6,
    get_listindex6,
    get_multiindex3,
    get_listindex3
export differentiate, poisson_bracket, changevars, evaluatepolynomial
export unit_variableset, terms_ofdegree

export find_gamma,
    get_hamiltonian_flow, eoms_polynomial, compose_polynomial_eoms, get_hamiltonian_jacobian

export c_n, dipole_expansion
export diagonalizing_matrix, complexifying_change, realify
export birkhoff_term, resonant_term, obtain_generating_function!, lie_series_transformation!
export normalform2birkhoff, birkhoff2normalform, normalform2resonant, resonant2normalform

export synodic2recentered,
    synodic2normalform,
    synodic2normalform_numerical,
    synodic2actionangle,
    synodic2actionangle_numerical,
    recentered2synodic,
    recentered2diagonal,
    diagonal2recentered,
    diagonal2normalform,
    diagonal2normalform_numerical,
    normalform2synodic,
    normalform2diagonal,
    normalform2actionangle,
    actionangle2synodic,
    actionangle2normalform,
    actionangle2pureactionangle,
    pureactionangle2actionangle



export generate_birkhoff_NF,
    generate_resonant_NF, generate_forwardchange, generate_inversechange

export compute_normalform, store_NF, compute_store_NF



end # module NFCR3BP
