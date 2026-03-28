function diagonalizing_matrix(sys::NamedTuple)
    lambda1 = sys.lambda1
    omega1 = sys.omega1
    omega2 = sys.omega2

    c2 = sys.dipole_coefficients[2]

    d1 = 2*lambda1*((4+3*c2)*lambda1^2+4+5*c2-6*c2^2)
    d2 = omega1*((4+3*c2)*omega1^2-4-5*c2+6*c2^2)
    s1 = sqrt(d1)
    s2 = sqrt(d2)

    C_tf = zeros(6, 6)

    # Row 1
    C_tf[1, 1] = 2*lambda1/s1
    C_tf[1, 4] = -2*lambda1/s1
    C_tf[1, 5] = 2*omega1/s2
    # Row 2
    C_tf[2, 1] = (lambda1^2-2*c2-1)/s1
    C_tf[2, 2] = (-omega1^2-2*c2-1)/s2
    C_tf[2, 4] = (lambda1^2-2*c2-1)/s1
    # Row 3
    C_tf[3, 3] = 1/sqrt(omega2)
    # Row 4
    C_tf[4, 1] = (lambda1^2+2*c2+1)/s1
    C_tf[4, 2] = (-omega1^2+2*c2+1)/s2
    C_tf[4, 4] = (lambda1^2+2*c2+1)/s1
    # Row 5
    C_tf[5, 1] = (lambda1^3+(1-2*c2)*lambda1)/s1
    C_tf[5, 4] = (-lambda1^3-(1-2*c2)*lambda1)/s1
    C_tf[5, 5] = (-omega1^3+(1-2*c2)*omega1)/s2
    # Row 6
    C_tf[6, 6] = sqrt(omega2)

    return C_tf
end

function diagonalizing_matrix(mu::Real, lagrangepoint::Integer)
    gamma = find_gamma(mu, lagrangepoint)
    sys = (mu = mu, gamma = gamma, lagrange_point = lagrangepoint)

    c2 = c_n(sys, 2)
    lambda1 = sqrt(0.5*(c2-2+sqrt(9*c2^2-8*c2)))
    omega1 = sqrt(-0.5*(c2-2-sqrt(9*c2^2-8*c2)))
    omega2 = sqrt(c2)

    d1 = 2*lambda1*((4+3*c2)*lambda1^2+4+5*c2-6*c2^2)
    d2 = omega1*((4+3*c2)*omega1^2-4-5*c2+6*c2^2)
    s1 = sqrt(d1)
    s2 = sqrt(d2)

    C_tf = zeros(6, 6)

    # Row 1
    C_tf[1, 1] = 2*lambda1/s1
    C_tf[1, 4] = -2*lambda1/s1
    C_tf[1, 5] = 2*omega1/s2
    # Row 2
    C_tf[2, 1] = (lambda1^2-2*c2-1)/s1
    C_tf[2, 2] = (-omega1^2-2*c2-1)/s2
    C_tf[2, 4] = (lambda1^2-2*c2-1)/s1
    # Row 3
    C_tf[3, 3] = 1/sqrt(omega2)
    # Row 4
    C_tf[4, 1] = (lambda1^2+2*c2+1)/s1
    C_tf[4, 2] = (-omega1^2+2*c2+1)/s2
    C_tf[4, 4] = (lambda1^2+2*c2+1)/s1
    # Row 5
    C_tf[5, 1] = (lambda1^3+(1-2*c2)*lambda1)/s1
    C_tf[5, 4] = (-lambda1^3-(1-2*c2)*lambda1)/s1
    C_tf[5, 5] = (-omega1^3+(1-2*c2)*omega1)/s2
    # Row 6
    C_tf[6, 6] = sqrt(omega2)

    return C_tf
end
