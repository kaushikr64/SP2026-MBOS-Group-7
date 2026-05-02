"""
    Abstract Polynomial type
"""
abstract type AbstractHomogenousPolynomial end
abstract type AbstractActionAnglePolynomial end
const POLYTOL = 1e-14

"""
    Homogenous Polynomials
"""

struct RealHomogenousPolynomial{T1<:Integer,T2<:SparseVector{<:Real}} <:
       AbstractHomogenousPolynomial
    order::T1
    terms::T2
end

struct ComplexHomogenousPolynomial{T1<:Integer,T2<:SparseVector{<:Complex}} <:
       AbstractHomogenousPolynomial
    order::T1
    terms::T2
end

struct BirkhoffActionAnglePolynomialDegree{
    T1<:Integer,
    T2<:SparseVector{<:Real}
}
    order::T1
    terms::T2
end


struct ResonantActionAnglePolynomialDegree{
    T1<:Integer,
    T2<:SparseVector{<:Real},
    T3<:SparseVector{<:SVector},
}
    order::T1
    integrable_terms::T2
    resonant_terms::T3
end


"""
    Polynomial
"""
struct MixedDegreePolynomial{T1<:AbstractVector{<:AbstractHomogenousPolynomial}}
    degrees::T1
end

struct ResonantActionAnglePolynomial{
    T1<:AbstractVector{<:ResonantActionAnglePolynomialDegree},
} <: AbstractActionAnglePolynomial
    degrees::T1
end

struct BirkhoffActionAnglePolynomial{
    T1<:AbstractVector{<:BirkhoffActionAnglePolynomialDegree},
} <: AbstractActionAnglePolynomial
    degrees::T1
end

"""
    Constructors
"""
function RealHomogenousPolynomial(coeff::Real, exponent::Vector{<:Integer})
    order = sum(exponent)
    num_monomials = get_num_monomials(6, order)
    listindex = get_listindex6(order, exponent)
    terms = sparsevec([listindex], [coeff], num_monomials)
    return RealHomogenousPolynomial(order, droptol!(terms, POLYTOL));
end

function ComplexHomogenousPolynomial(coeff::Complex, exponent::Vector{<:Integer})
    order = sum(exponent)
    num_monomials = get_num_monomials(6, order)
    listindex = get_listindex6(order, exponent)
    terms = sparsevec([listindex], [coeff], num_monomials)
    return ComplexHomogenousPolynomial(order, droptol!(terms, POLYTOL));
end


function RealHomogenousPolynomial(
    coeffs::Vector{<:Real},
    exponents::Vector{<:Vector{<:Integer}},
)
    order = sum(exponents[1])
    num_monomials = get_num_monomials(6, order)
    listindices = Vector{Integer}(undef, length(order))
    for i in eachindex(listindices)
        listindices[i] = get_listindex6(order, exponents[i])
    end
    terms = sparsevec(listindices, coeffs, num_monomials)
    return RealHomogenousPolynomial(order, droptol!(terms, POLYTOL));
end

function ComplexHomogenousPolynomial(
    coeffs::Vector{<:Complex},
    exponents::Vector{<:Vector{<:Integer}},
)
    order = sum(exponents[1])
    num_monomials = get_num_monomials(6, order)
    listindices = Vector{Integer}(undef, length(order))
    for i in eachindex(listindices)
        listindices[i] = get_listindex6(order, exponents[i])
    end
    terms = sparsevec(listindices, coeffs, num_monomials)
    return ComplexHomogenousPolynomial(order, droptol!(terms, POLYTOL));
end

function ZeroHomogenousPolynomial(order::Integer)
    return RealHomogenousPolynomial(Float64(0), [order, 0, 0, 0, 0, 0])
end


function CombinePolynomials(constituents::Vector{<:AbstractHomogenousPolynomial})
    max_polynomialorder = MAXDEGREE
    polynomial_vector = Vector{AbstractHomogenousPolynomial}(undef, max_polynomialorder+1)
    for i = 0:max_polynomialorder
        polynomial_vector[i+1] = ZeroHomogenousPolynomial(i)
        for ii in eachindex(constituents)
            if constituents[ii].order == i
                polynomial_vector[i+1] += constituents[ii]
            end
        end
    end
    return MixedDegreePolynomial(polynomial_vector)
end

