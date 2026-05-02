function generate_birkhoff_NF(mu::Real, lagrangepoint::Integer, maxorder::Integer)
    # Compute gamma from the given value of μ
    gamma = find_gamma(mu, lagrangepoint)

    #Term removal criteria
    remove_terms(multiindex) = !birkhoff_term(multiindex)

    # Initialize the system parameters dictionary
    sys = (
        mu = mu,
        gamma = gamma,
        maxorder_NF = maxorder,
        removal_criteria = remove_terms,
        lagrange_point = lagrangepoint,
    )

    # Compute dipole expansion coefficients
    cn = c_n(sys)
    c2 = cn[2]

    # Compute system eigenvalues (normal modes)
    lambda1 = sqrt(0.5*(c2-2+sqrt(9*c2^2-8*c2)))
    omega1 = sqrt(-0.5*(c2-2-sqrt(9*c2^2-8*c2)))
    omega2 = sqrt(c2)
    lambda = [lambda1; 1im*omega1; 1im*omega2]

    # Push frequencies to system parameters dictionary
    sys = merge(
        sys,
        (
            dipole_coefficients = cn,
            lambda1 = lambda1,
            omega1 = omega1,
            omega2 = omega2,
            eigenvalues = lambda,
        ),
    )

    # Compute diagonalizing change, obtaining the recentered variables as a function of 
    # the complex diagonalized variables
    recentered_vars = diagonalizing_matrix(sys)*complexifying_change()*unit_variableset()
    xbar = recentered_vars[1]
    ybar = recentered_vars[2]
    zbar = recentered_vars[3]
    pxbar = recentered_vars[4]
    pybar = recentered_vars[5]
    pzbar = recentered_vars[6]

    # Now obtain the complex diagonal Hamiltonian
    non_dipole_terms = 0.5*(pxbar^2+pybar^2+pzbar^2)-xbar*pybar+ybar*pxbar
    dipole_terms = dipole_expansion(recentered_vars, sys)
    H_complex_diagonal = non_dipole_terms - dipole_terms

    # Run a copy of this Hamiltonian through the normalization process
    H_complex_normal = deepcopy(H_complex_diagonal)
    H_complex_normal, G_complex = lie_series_transformation!(H_complex_normal, sys)

    # Realify both H and G
    G_real = realify(G_complex)
    H_real_normal = realify(H_complex_normal)


    # Finally, compute the action angle Hamiltonian
    H_action_angle = normalform2birkhoff(H_complex_normal)

    return H_complex_normal, H_real_normal, G_real, H_action_angle
end


function generate_resonant_NF(mu::Real, lagrangepoint::Integer, maxorder::Integer)
    # Compute gamma from the given value of μ
    gamma = find_gamma(mu, lagrangepoint)

    #Term removal criteria
    remove_terms(multiindex) = !resonant_term(multiindex)

    # Initialize the system parameters dictionary
    sys = (
        mu = mu,
        gamma = gamma,
        maxorder_NF = maxorder,
        removal_criteria = remove_terms,
        lagrange_point = lagrangepoint,
    )

    # Compute dipole expansion coefficients
    cn = c_n(sys)
    c2 = cn[2]

    # Compute system eigenvalues (normal modes)
    lambda1 = sqrt(0.5*(c2-2+sqrt(9*c2^2-8*c2)))
    omega1 = sqrt(-0.5*(c2-2-sqrt(9*c2^2-8*c2)))
    omega2 = sqrt(c2)
    lambda = [lambda1; 1im*omega1; 1im*omega2]

    # Push frequencies to system parameters dictionary
    sys = merge(
        sys,
        (
            dipole_coefficients = cn,
            lambda1 = lambda1,
            omega1 = omega1,
            omega2 = omega2,
            eigenvalues = lambda,
        ),
    )

    # Compute diagonalizing change, obtaining the recentered variables as a function of 
    # the complex diagonalized variables
    recentered_vars = diagonalizing_matrix(sys)*complexifying_change()*unit_variableset()
    xbar = recentered_vars[1]
    ybar = recentered_vars[2]
    zbar = recentered_vars[3]
    pxbar = recentered_vars[4]
    pybar = recentered_vars[5]
    pzbar = recentered_vars[6]

    # Now obtain the complex diagonal Hamiltonian
    non_dipole_terms = 0.5*(pxbar^2+pybar^2+pzbar^2)-xbar*pybar+ybar*pxbar
    dipole_terms = dipole_expansion(recentered_vars, sys)
    H_complex_diagonal = non_dipole_terms - dipole_terms

    # Run a copy of this Hamiltonian through the normalization process
    H_complex_normal = deepcopy(H_complex_diagonal)
    H_complex_normal, G_complex = lie_series_transformation!(H_complex_normal, sys)

    # Realify both H and G
    G_real = realify(G_complex)
    H_real_normal = realify(H_complex_normal)


    # Finally, compute the action angle Hamiltonian
    H_action_angle = normalform2resonant(H_complex_normal)

    return H_complex_normal, H_real_normal, G_real, H_action_angle
end


function generate_RCM(mu::Real, lagrangepoint::Integer, maxorder::Integer)
    # Compute gamma from the given value of μ
    gamma = find_gamma(mu, lagrangepoint)

    #Term removal criteria
    remove_terms(multiindex) = rcm_term(multiindex)

    # Initialize the system parameters dictionary
    sys = (
        mu = mu,
        gamma = gamma,
        maxorder_NF = maxorder,
        removal_criteria = remove_terms,
        lagrange_point = lagrangepoint,
    )

    # Compute dipole expansion coefficients
    cn = c_n(sys)
    c2 = cn[2]

    # Compute system eigenvalues (normal modes)
    lambda1 = sqrt(0.5*(c2-2+sqrt(9*c2^2-8*c2)))
    omega1 = sqrt(-0.5*(c2-2-sqrt(9*c2^2-8*c2)))
    omega2 = sqrt(c2)
    lambda = [lambda1; 1im*omega1; 1im*omega2]

    # Push frequencies to system parameters dictionary
    sys = merge(
        sys,
        (
            dipole_coefficients = cn,
            lambda1 = lambda1,
            omega1 = omega1,
            omega2 = omega2,
            eigenvalues = lambda,
        ),
    )

    # Compute diagonalizing change, obtaining the recentered variables as a function of 
    # the complex diagonalized variables
    recentered_vars = diagonalizing_matrix(sys)*complexifying_change()*unit_variableset()
    xbar = recentered_vars[1]
    ybar = recentered_vars[2]
    zbar = recentered_vars[3]
    pxbar = recentered_vars[4]
    pybar = recentered_vars[5]
    pzbar = recentered_vars[6]

    # Now obtain the complex diagonal Hamiltonian
    non_dipole_terms = 0.5*(pxbar^2+pybar^2+pzbar^2)-xbar*pybar+ybar*pxbar
    dipole_terms = dipole_expansion(recentered_vars, sys)
    H_complex_diagonal = non_dipole_terms - dipole_terms

    # Run a copy of this Hamiltonian through the normalization process
    H_complex_normal = deepcopy(H_complex_diagonal)
    H_complex_normal, G_complex = lie_series_transformation!(H_complex_normal, sys)

    # Realify both H and G
    G_real = realify(G_complex)
    H_real_normal = realify(H_complex_normal)

    return H_complex_normal, H_real_normal, G_real
end


function generate_forwardchange(G, maxorder)
    variable_transformations = Vector{MixedDegreePolynomial}(undef, 6)
    variableset = unit_variableset()
    Threads.@threads for i = 1:6
        variable_transformation = CombinePolynomials([variableset[i]]) # Initialize it as multi degree
        for order = 3:maxorder # The orders of the generating function

            order_increment = order-2
            generating_function = terms_ofdegree(G, order)
            for ii = maxorder:-1:(order-1)
                degreeindex = ii+1
                stepstaken = maxorder-ii
                substepnum = div(stepstaken, order_increment) + 1
                poiss_brac = terms_ofdegree(variable_transformation, ii-order_increment)
                for ii = 1:substepnum
                    poiss_brac = poisson_bracket(poiss_brac, generating_function)
                    variable_transformation.degrees[degreeindex+(ii-1)*order_increment] +=
                        ((1/factorial(ii))*poiss_brac)
                end
            end
        end
        variable_transformations[i] = variable_transformation
    end
    return variable_transformations
end

function generate_inversechange(G, maxorder)
    variable_transformations = Vector{MixedDegreePolynomial}(undef, 6)
    variableset = unit_variableset()
    Threads.@threads for i = 1:6
        variable_transformation = CombinePolynomials([variableset[i]]) # Initialize it as multi degree
        for order = maxorder:-1:3 # The orders of the generating function
            order_increment = order-2
            generating_function = -1*terms_ofdegree(G, order)
            for ii = maxorder:-1:(order-1)
                degreeindex = ii+1
                stepstaken = maxorder-ii
                substepnum = div(stepstaken, order_increment) + 1
                poiss_brac = terms_ofdegree(variable_transformation, ii-order_increment)
                for ii = 1:substepnum
                    poiss_brac = poisson_bracket(poiss_brac, generating_function)
                    variable_transformation.degrees[degreeindex+(ii-1)*order_increment] +=
                        ((1/factorial(ii))*poiss_brac)
                end
            end
        end
        variable_transformations[i] = variable_transformation
    end
    return variable_transformations
end

