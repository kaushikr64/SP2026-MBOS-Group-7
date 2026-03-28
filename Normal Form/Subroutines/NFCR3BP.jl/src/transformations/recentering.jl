function c_n(sys::NamedTuple, n::Int64)
    mu = sys.mu
    gamma = sys.gamma
    lagrangepoint = sys.lagrange_point
    if lagrangepoint == 1
        c = (1/gamma^3) * (mu + ((-1)^n) * ((1-mu) * (gamma^(n+1))) / (1-gamma)^(n+1))
    elseif lagrangepoint == 2
        c =
            (1/gamma^3) *
            (((-1)^n)*mu + ((-1)^n) * ((1-mu) * (gamma^(n+1))) / (1+gamma)^(n+1))
    elseif lagrangepoint == 3
        c = (((-1)^n)/gamma^3)*(1-mu+(mu*gamma^(n+1))/(1+gamma)^(n+1))
    end

    return c
end

function c_n(sys::NamedTuple)
    mu = sys.mu
    gamma = sys.gamma
    lagrangepoint = sys.lagrange_point
    n = 1:sys.maxorder_NF
    if lagrangepoint == 1
        c = [
            (1/gamma^3) * (mu + ((-1)^i) * ((1-mu) * (gamma^(i+1))) / (1-gamma)^(i+1))
            for i in n
        ]
    elseif lagrangepoint == 2
        c = [
            (1/gamma^3) *
            (((-1)^i)*mu + ((-1)^i) * ((1-mu) * (gamma^(i+1))) / (1+gamma)^(i+1)) for
            i in n
        ]
    elseif lagrangepoint == 3
        c = [(((-1)^i)/gamma^3)*(1-mu+(mu*gamma^(i+1))/(1+gamma)^(i+1)) for i in n]
    end

    return c
end

function dipole_expansion(state_vars::Vector{<:AbstractHomogenousPolynomial}, sys)
    mu = sys.mu
    gamma = sys.gamma
    n = sys.maxorder_NF
    cn = sys.dipole_coefficients

    # variables
    x = state_vars[1]
    y = state_vars[2]
    z = state_vars[3]

    # Obtain recurrent polynomials
    T = Vector{Any}(undef, n+1)
    T[1] = 1
    T[2] = x
    for i = 2:n
        T[i+1] = ((2*i-1)/i) * x * T[i] - ((i-1)/i) * (x^2 + y^2 + z^2) * T[i-1]
    end

    A = cn[2:end] .* T[3:end]

    return CombinePolynomials(A)
end

function recenter_rescaling_matrices(mu::Real, lagrange_point::Integer)
    gamma = find_gamma(mu, lagrange_point)

    if lagrange_point < 3
        a = 1+(-1)^lagrange_point*gamma
        delta = 1
    elseif lagrange_point == 3
        a = -gamma
        delta = -1
    end
    recentering = SVector(a-mu, 0, 0, 0, a-mu, 0)

    # Rescaling
    rescaling = Diagonal([delta*gamma, delta*gamma, gamma, delta*gamma, delta*gamma, gamma])

    # Convert q,qdot to q,p    
    qqdot2qp = Matrix(1I(6))
    qqdot2qp[4, 2] = 1
    qqdot2qp[5, 1] = -1

    return rescaling, qqdot2qp, recentering

end
