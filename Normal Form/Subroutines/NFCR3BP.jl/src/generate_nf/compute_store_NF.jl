"""
    Compute and Store normal form data of desired type to desired path
"""
function compute_normalform(
    mu::Real,
    lagrangepoint::Integer,
    normalization_order::Integer,
    normalization_type::String,
)
    # Obtain the common linear changes
    rescaling, qqdot2qp, recentering = recenter_rescaling_matrices(mu, lagrangepoint)
    diagonalizing = diagonalizing_matrix(mu, lagrangepoint)

    if normalization_type == "Birkhoff"
        println("Normalizing Hamiltonian up to order ", normalization_order)
        time_NF = @elapsed (H_complex_normal, H_real_normal, G_real, H_action_angle) =
            generate_birkhoff_NF(mu, lagrangepoint, normalization_order)
        println("Birkhoff normal form computed in ", time_NF, " seconds.")

        nf2aa = normalform2birkhoff
        aa2nf = birkhoff2normalform

        hamiltonians = (
            complex_normal = H_complex_normal,
            real_normal = H_real_normal,
            action_angle = H_action_angle,
        )
    elseif normalization_type == "Resonant"
        println("Normalizing Hamiltonian up to order ", normalization_order)
        time_NF = @elapsed (H_complex_normal, H_real_normal, G_real, H_action_angle) =
            generate_resonant_NF(mu, lagrangepoint, normalization_order)
        println("Resonant normal form computed in ", time_NF, " seconds.")

        nf2aa = normalform2resonant
        aa2nf = resonant2normalform

        hamiltonians = (
            complex_normal = H_complex_normal,
            real_normal = H_real_normal,
            action_angle = H_action_angle,
        )
    end

    time_eoms = @elapsed normalform_flow = get_hamiltonian_flow(H_real_normal)
    println("Equations of motion computed in ", time_eoms, " seconds.")

    time_jac = @elapsed normalform_jac = get_hamiltonian_jacobian(normalform_flow)
    println("Jacobian computed in ", time_eoms, " seconds.")

    println("Computing forward change")
    time_forwards = @elapsed forward_transformation =
        generate_forwardchange(G_real, normalization_order)
    println("Forward change computed in ", time_forwards, " seconds.")

    time_backwards = @elapsed backward_transformation =
        generate_inversechange(G_real, normalization_order)
    println("Backward change computed in ", time_backwards, " seconds.")

    time_G_flow = @elapsed G_flow = get_hamiltonian_flow(G_real)
    println("Generating function flow computed in ", time_G_flow, " seconds.")

    time_H_flow = @elapsed H_flow = get_hamiltonian_flow(H_action_angle)
    println("Hamiltonian flow computed in ", time_H_flow, " seconds.")


    generating_function = G_real
    derivatives = (
        NF_flow = normalform_flow::Vector{MixedDegreePolynomial},
        jac = normalform_jac::Matrix{MixedDegreePolynomial},
        pureAA_flow = H_flow::Vector{ResonantActionAnglePolynomial},
        AA_flow = Vector{Union{MixedDegreePolynomial,ResonantActionAnglePolynomial}}(
            vcat(normalform_flow[1:2], H_flow[3:6]),
        ),
    )
    transformations = (
        recentering = recentering,
        qqdot2qp = qqdot2qp,
        rescaling = rescaling,
        diagonalizing = diagonalizing,
        forward = forward_transformation,
        backward = backward_transformation,
        G_flow = G_flow,
        nf2aa = nf2aa,
        aa2nf = aa2nf,
    )

    return hamiltonians, generating_function, derivatives, transformations
end

function store_NF(
    hamiltonians,
    generating_function,
    derivatives,
    transformations,
    storage_file::String,
)
    @save storage_file hamiltonians generating_function derivatives transformations
    return nothing
end

function compute_store_NF(storage_file::String, params::NamedTuple)
    hamiltonians, generating_function, derivatives, transformations = compute_normalform(
        params.mu,
        params.lagrange_point,
        params.normalization_order,
        params.normalization_type,
    )

    store_NF(hamiltonians, generating_function, derivatives, transformations, storage_file)
    return nothing
end
