function propagate(f::Function, x0, tf::Real, p; method = DP8(), kwargs...)
    prob = ODEProblem(f, x0, (0, tf), p; kwargs...)
    sol = solve(prob, method)
    return sol(tf)
end

function eoms_polynomial(xi, polynomial_dxidt, t)
    dxidt = Vector{eltype(xi)}(undef, length(xi))
    Threads.@threads for i = 1:length(xi)
        dxidt[i] = evaluatepolynomial(polynomial_dxidt[i], xi)
    end
    return SVector{length(xi)}(dxidt)
end

function eoms_polynomial(
    xi_aa,
    polynomial_dxidt::Vector{Union{MixedDegreePolynomial,ResonantActionAnglePolynomial}},
    t,
)

    n = length(polynomial_dxidt)
    xi_nf = resonant2normalform(xi_aa)
    dxidt = Vector{eltype(xi_aa)}(undef, n)

    Threads.@threads for i = 1:n
        p = polynomial_dxidt[i]
        if p isa MixedDegreePolynomial
            dxidt[i] = evaluatepolynomial(p, xi_nf)
        else
            dxidt[i] = evaluatepolynomial(p, xi_aa)
        end
    end
    return SVector{n}(dxidt)
end

function compose_polynomial_eoms(
    vector_flow::Vector{MixedDegreePolynomial},
    matrix_flow::Matrix{MixedDegreePolynomial},
)
    return vcat(vector_flow, vec(matrix_flow))
end
