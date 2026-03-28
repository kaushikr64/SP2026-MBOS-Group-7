function L1_quintic(mu::Real)
    f(γ) = γ^5 - (3 - mu) * γ^4 + (3 - 2 * mu) * γ^3 - mu * γ^2 + 2 * mu * γ - mu
    f_dash(γ) = 5 * γ^4 - 4 * (3 - mu) * γ^3 + 3 * (3 - 2 * mu) * γ^2 - 2 * mu * γ + 2 * mu
    return f, f_dash
end

function L2_quintic(mu::Real)
    f(γ) = γ^5 + (3 - mu) * γ^4 + (3 - 2 * mu) * γ^3 - mu * γ^2 - 2 * mu * γ - mu
    f_dash(γ) = 5 * γ^4 + 4 * (3 - mu) * γ^3 + 3 * (3 - 2 * mu) * γ^2 - 2 * mu * γ - 2 * mu
    return f, f_dash
end

function L3_quintic(mu::Real)
    f(γ) =
        γ^5 + (2 + mu) * γ^4 + (1 + 2 * mu) * γ^3 - (1 - mu) * γ^2 - 2 * (1 - mu) * γ -
        (1 - mu);
    f_dash(γ) =
        5 * γ^4 + 4 * (2 + mu) * γ^3 + 3 * (1 + 2 * mu) * γ^2 - 2 * (1 - mu) * γ -
        2 * (1 - mu);
    return f, f_dash
end

function newton_raphson(γ_0::Real, f::Function, f_dash::Function; solver_tol = 1e-16::Real)
    converged = false
    γ_n = γ_0
    while !converged
        γ_n1 = γ_n - f(γ_n)/f_dash(γ_n)
        if abs(γ_n1-γ_n) < solver_tol
            return γ_n
        else
            γ_n = γ_n1
        end
    end
end

function find_gamma(mu::Real, lagrangepoint::Integer; solver_tol = 1e-16::Real)
    if lagrangepoint == 1
        f, f_dash = L1_quintic(mu)
    elseif lagrangepoint == 2
        f, f_dash = L2_quintic(mu)
    elseif lagrangepoint == 3
        f, f_dash = L3_quintic(mu)
    end
    return newton_raphson(0, f, f_dash; solver_tol = solver_tol)
end
