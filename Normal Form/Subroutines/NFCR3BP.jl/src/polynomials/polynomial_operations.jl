import Base: +, -, *, ^, zero

"""
    Define zero
"""
function zero(a::AbstractHomogenousPolynomial)
    return RealHomogenousPolynomial(0, [a.order, 0, 0, 0, 0, 0])
end

"""
    Addition & subtraction
"""
function +(a::RealHomogenousPolynomial, b::RealHomogenousPolynomial)
    c = RealHomogenousPolynomial(a.order, droptol!(a.terms + b.terms, POLYTOL))
    return c
end

function +(a::AbstractHomogenousPolynomial, b::AbstractHomogenousPolynomial)
    c = ComplexHomogenousPolynomial(
        a.order,
        SparseVector{Complex}(droptol!(a.terms+b.terms, POLYTOL)),
    )
    return c
end

function -(a::RealHomogenousPolynomial, b::RealHomogenousPolynomial)
    c = RealHomogenousPolynomial(a.order, droptol!(a.terms - b.terms, POLYTOL))
    return c
end

function -(a::AbstractHomogenousPolynomial, b::AbstractHomogenousPolynomial)
    c = ComplexHomogenousPolynomial(
        a.order,
        SparseVector{Complex}(droptol!(a.terms-b.terms, POLYTOL)),
    )
    return c
end

function +(a::MixedDegreePolynomial, b::MixedDegreePolynomial)
    c = Vector{AbstractHomogenousPolynomial}(undef, MAXDEGREE)
    for i = 1:(MAXDEGREE+1)
        c[i] = a.degrees[i]+b.degrees[i]
    end
    return MixedDegreePolynomial(c)
end

function +(a::AbstractHomogenousPolynomial, b::MixedDegreePolynomial)
    c_degrees = copy(b.degrees)
    c_degrees[a.order+1] = a + b.degrees[a.order+1]
    return CombinePolynomials(c_degrees)
end

function +(a::MixedDegreePolynomial, b::AbstractHomogenousPolynomial)
    c_degrees = copy(b.degrees)
    c_degrees[b.order+1] = a.degrees[b.order+1] + b
    return CombinePolynomials(c_degrees)
end

function +(a::ResonantActionAnglePolynomialDegree, b::ResonantActionAnglePolynomialDegree)
    c_integrable_terms = a.integrable_terms + b.integrable_terms
    a_nonzero = rowvals(a.resonant_terms)
    b_nonzero = rowvals(b.resonant_terms)
    all_indices = sort(union(a_nonzero, b_nonzero))
    c_resonant_idxs = Int[]
    c_resonant_vals = eltype(a.resonant_terms)[]
    for i in all_indices
        if i in a_nonzero && i in b_nonzero
            a_v = a.resonant_terms[i]
            b_v = b.resonant_terms[i]
            new_num = vcat(a_v[1], b_v[1])
            new_trig1 = vcat(a_v[2], b_v[2])
            new_trig2 = vcat(a_v[3], b_v[3])
            keep = abs.(new_num) .> POLYTOL
            new_num, new_trig1, new_trig2 = new_num[keep], new_trig1[keep], new_trig2[keep]
            if !isempty(new_num)
                push!(c_resonant_idxs, i)
                push!(c_resonant_vals, SVector(new_num, new_trig1, new_trig2))
            end
        elseif i in a_nonzero
            push!(c_resonant_idxs, i)
            push!(c_resonant_vals, a.resonant_terms[i])
        else
            push!(c_resonant_idxs, i)
            push!(c_resonant_vals, b.resonant_terms[i])
        end
    end
    c_resonant_terms = sparsevec(c_resonant_idxs, c_resonant_vals, length(a.resonant_terms))
    return ResonantActionAnglePolynomialDegree(
        a.order,
        droptol!(c_integrable_terms, POLYTOL),
        c_resonant_terms,
    )
end

function +(a::ResonantActionAnglePolynomial, b::ResonantActionAnglePolynomial)
    c = Vector{ResonantActionAnglePolynomialDegree}(undef, length(a.degrees))
    for i = 1:length(a.degrees)
        c[i] = a.degrees[i]+b.degrees[i]
    end
    return ResonantActionAnglePolynomial(c)
end

function -(a::MixedDegreePolynomial, b::MixedDegreePolynomial)
    c = Vector{AbstractHomogenousPolynomial}(undef, MAXDEGREE)
    for i = 1:(MAXDEGREE+1)
        c[i] = a.degrees[i]-b.degrees[i]
    end
    return MixedDegreePolynomial(c)
end

function -(a::AbstractHomogenousPolynomial, b::MixedDegreePolynomial)
    c_degrees = -1*b.degrees
    c_degrees[a.order+1] = a - b.degrees[a.order+1]
    return CombinePolynomials(c_degrees)
end

function -(a::MixedDegreePolynomial, b::AbstractHomogenousPolynomial)
    c_degrees = -1*a.degrees
    c_degrees[b.order+1] = a.degrees[b.order+1] - b
    return CombinePolynomials(c_degrees)
end

function -(a::ResonantActionAnglePolynomialDegree, b::ResonantActionAnglePolynomialDegree)
    c_integrable_terms = a.integrable_terms - b.integrable_terms
    a_nonzero = rowvals(a.resonant_terms)
    b_nonzero = rowvals(b.resonant_terms)
    all_indices = sort(union(a_nonzero, b_nonzero))
    c_resonant_idxs = Int[]
    c_resonant_vals = eltype(a.resonant_terms)[]
    for i in all_indices
        if i in a_nonzero && i in b_nonzero
            a_v = a.resonant_terms[i]
            b_v = b.resonant_terms[i]
            new_num = vcat(a_v[1], -b_v[1])
            new_trig1 = vcat(a_v[2], b_v[2])
            new_trig2 = vcat(a_v[3], b_v[3])
            keep = abs.(new_num) .> POLYTOL
            new_num, new_trig1, new_trig2 = new_num[keep], new_trig1[keep], new_trig2[keep]
            if !isempty(new_num)
                push!(c_resonant_idxs, i)
                push!(c_resonant_vals, SVector(new_num, new_trig1, new_trig2))
            end
        elseif i in a_nonzero
            push!(c_resonant_idxs, i)
            push!(c_resonant_vals, a.resonant_terms[i])
        else
            b_v = b.resonant_terms[i]
            push!(c_resonant_idxs, i)
            push!(c_resonant_vals, SVector(-b_v[1], b_v[2], b_v[3]))
        end
    end
    c_resonant_terms = sparsevec(c_resonant_idxs, c_resonant_vals, length(a.resonant_terms))
    return ResonantActionAnglePolynomialDegree(
        a.order,
        droptol!(c_integrable_terms, POLYTOL),
        c_resonant_terms,
    )
end

function -(a::ResonantActionAnglePolynomial, b::ResonantActionAnglePolynomial)
    c = Vector{ResonantActionAnglePolynomialDegree}(undef, length(a.degrees))
    for i = 1:length(a.degrees)
        c[i] = a.degrees[i]-b.degrees[i]
    end
    return ResonantActionAnglePolynomial(c)
end

"""
    Multiplication
"""
function *(a::RealHomogenousPolynomial, b::RealHomogenousPolynomial)
    # set up vectors
    c_order = a.order+b.order
    num_monomials = get_num_monomials(6, c_order)
    c_indices = Integer[]
    c_coeffs = Real[]

    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)
    b_indices = rowvals(b.terms)
    b_coeffs = nonzeros(b.terms)


    for i in eachindex(a_indices)
        for ii in eachindex(b_indices)
            a_index = a_indices[i]
            b_index = b_indices[ii]
            a_coeff = a_coeffs[i]
            b_coeff = b_coeffs[ii]
            a_multiindex = get_multiindex6(a.order, a_index)
            b_multiindex = get_multiindex6(b.order, b_index)

            # Actually do the multiplication
            c_coeff = a_coeff*b_coeff
            c_multiindex = a_multiindex+b_multiindex
            if sum(c_multiindex) <= MAXDEGREE
                c_index = get_listindex6(c_order, c_multiindex)

                push!(c_coeffs, c_coeff)
                push!(c_indices, c_index)
            end
        end
    end

    c_terms = sparsevec(c_indices, c_coeffs, num_monomials)

    return RealHomogenousPolynomial(c_order, c_terms)
end

function *(a::AbstractHomogenousPolynomial, b::AbstractHomogenousPolynomial)
    # set up vectors
    c_order = a.order+b.order
    num_monomials = get_num_monomials(6, c_order)
    c_indices = Integer[]
    c_coeffs = Complex[]

    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)
    b_indices = rowvals(b.terms)
    b_coeffs = nonzeros(b.terms)


    for i in eachindex(a_indices)
        for ii in eachindex(b_indices)
            a_index = a_indices[i]
            b_index = b_indices[ii]
            a_coeff = a_coeffs[i]
            b_coeff = b_coeffs[ii]
            a_multiindex = get_multiindex6(a.order, a_index)
            b_multiindex = get_multiindex6(b.order, b_index)

            # Actually do the multiplication
            c_coeff = a_coeff*b_coeff
            c_multiindex = a_multiindex+b_multiindex
            if sum(c_multiindex) <= MAXDEGREE
                c_index = get_listindex6(c_order, c_multiindex)

                push!(c_coeffs, c_coeff)
                push!(c_indices, c_index)
            end
        end
    end

    c_terms = sparsevec(c_indices, c_coeffs, num_monomials)

    return ComplexHomogenousPolynomial(c_order, c_terms)
end

function *(a::Real, b::RealHomogenousPolynomial)
    # set up vectors
    c_terms = SparseVector{Real}(a*b.terms)
    return RealHomogenousPolynomial(b.order, c_terms)
end

function *(a::RealHomogenousPolynomial, b::Real)
    return b*a
end

function *(a::Number, b::AbstractHomogenousPolynomial)
    # set up vectors
    c_terms = SparseVector{Complex}(a*b.terms)
    return ComplexHomogenousPolynomial(b.order, c_terms)
end

function *(a::AbstractHomogenousPolynomial, b::Number)
    return b*a
end

function *(a::Matrix, b::Vector{<:AbstractHomogenousPolynomial})
    return [sum(a[i, :] .* b) for i = 1:length(b)]
end

function *(a::Real, b::MixedDegreePolynomial)
    return MixedDegreePolynomial([a*polynomial for polynomial in b.degrees])
end

function *(a::Real, b::ResonantActionAnglePolynomial)
    return ResonantActionAnglePolynomial([a*polynomial for polynomial in b.degrees])
end

function *(a::Real, b::ResonantActionAnglePolynomialDegree)
    return ResonantActionAnglePolynomialDegree(
        b.order,
        sparsevec(
            rowvals(b.integrable_terms),
            Vector{Real}([a*coefficient for coefficient in nonzeros(b.integrable_terms)]),
            length(b.integrable_terms),
        ),
        sparsevec(
            rowvals(b.resonant_terms),
            Vector{SVector{3}}([
                SVector(a.*coefficient[1], coefficient[2], coefficient[3]) for
                coefficient in nonzeros(b.resonant_terms)
            ]),
            length(b.resonant_terms),
        ),
    )
end


"""
    Power
"""
function ^(a::AbstractHomogenousPolynomial, b::Integer)
    exponentiated_value = 1
    for i = 1:b
        exponentiated_value *= a
    end
    return (exponentiated_value)
end


"""
    Derivative
"""
function differentiate(a::RealHomogenousPolynomial, b::Integer)
    if a.order == 0
        return ZeroHomogenousPolynomial(0)
    end

    output_order = a.order-1

    output_terms = spzeros(get_num_monomials(6, output_order))

    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)


    for i in eachindex(a_indices)
        exponent = get_multiindex6(a.order, a_indices[i])
        if exponent[b]!=0
            output_coeff = a_coeffs[i]*exponent[b]
            output_exponent = Vector(exponent)
            output_exponent[b] -= 1
            output_listindex = get_listindex6(output_order, output_exponent)
            output_terms[output_listindex] = output_coeff
        end
    end
    return RealHomogenousPolynomial(output_order, droptol!(output_terms, POLYTOL))
end

function differentiate(a::ComplexHomogenousPolynomial, b::Integer)
    output_order = a.order-1
    output_terms = spzeros(Complex, get_num_monomials(6, output_order))

    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)


    for i in eachindex(a_indices)
        exponent = get_multiindex6(a.order, a_indices[i])
        if exponent[b]!=0
            output_coeff = a_coeffs[i]*exponent[b]
            output_exponent = Vector(exponent)
            output_exponent[b] -= 1
            output_listindex = get_listindex6(output_order, output_exponent)
            output_terms[output_listindex] = output_coeff
        end
    end
    return ComplexHomogenousPolynomial(output_order, droptol!(output_terms, POLYTOL))
end

function differentiate(a::MixedDegreePolynomial, b::Integer)
    output_degrees = AbstractHomogenousPolynomial[]
    for degree in a.degrees
        if nnz(degree.terms) != 0
            push!(output_degrees, differentiate(degree, b))
        end
    end
    return CombinePolynomials(output_degrees)
end

function differentiate(a::ResonantActionAnglePolynomialDegree, b::Integer)
    a_integrable_indices = rowvals(a.integrable_terms)
    a_resonant_indices = rowvals(a.resonant_terms)
    a_integrable_coeffs = nonzeros(a.integrable_terms)
    a_resonant_coeffs = nonzeros(a.resonant_terms)
    if b < 4
        output_order = a.order-1
        output_integrable_terms = spzeros(Real, get_num_monomials(3, output_order))
        output_resonant_terms = spzeros(SVector{3}, get_num_monomials(3, output_order))
        # Derivative wrt action terms
        for i in eachindex(a_integrable_indices)
            exponent = get_multiindex3(a.order, a_integrable_indices[i])
            if exponent[b]!=0
                a_coeff = a_integrable_coeffs[i]
                output_coeff = a_coeff*exponent[b]
                output_exponent = Vector(exponent)
                output_exponent[b] -= 1
                output_listindex = get_listindex3(output_order, output_exponent)
                output_integrable_terms[output_listindex] = output_coeff
            end
        end

        for i in eachindex(a_resonant_indices)
            exponent = get_multiindex3(a.order, a_resonant_indices[i])
            if exponent[b]!=0
                a_coeff = a_resonant_coeffs[i]
                output_coeff = a_coeff[1]*exponent[b]
                output_exponent = Vector(exponent)
                output_exponent[b] -= 1
                output_listindex = get_listindex3(output_order, output_exponent)
                output_resonant_terms[output_listindex] =
                    SVector(output_coeff, a_coeff[2], a_coeff[3])
            end
        end
    else
        output_order = a.order
        output_integrable_terms = spzeros(Real, get_num_monomials(3, output_order))
        output_resonant_terms = spzeros(SVector{3}, get_num_monomials(3, output_order))
        # Derivative wrt resonant angle difference
        for i in eachindex(a_resonant_indices)
            listindex = a_resonant_indices[i]
            a_coeff = a_resonant_coeffs[i]
            a_numcoeff = a_coeff[1]
            a_sincoeff = a_coeff[2]
            a_coscoeff = a_coeff[3]
            if a_coscoeff != 0
                output_coeff = a_numcoeff.*a_coscoeff
                output_resonant_terms[listindex] = SVector(-output_coeff, a_coscoeff, a_sincoeff)
            elseif a_sincoeff != 0
                output_coeff = a_numcoeff.*a_sincoeff
                output_resonant_terms[listindex] = SVector(output_coeff, a_coscoeff, a_sincoeff)
            end

        end
    end
    return ResonantActionAnglePolynomialDegree(
        output_order,
        droptol!(output_integrable_terms, POLYTOL),
        output_resonant_terms,
    )
end

function differentiate(a::ResonantActionAnglePolynomial, b::Integer)
    if b < 4
        # Action derivative reduces order by 1; iterate from order-1 onwards (skip order-0 constant)
        output_degrees =
            Vector{ResonantActionAnglePolynomialDegree}(undef, length(a.degrees) - 1)
        for (out_idx, degree) in enumerate(a.degrees[2:end])
            output_degrees[out_idx] = differentiate(degree, b)
        end
    else
        # Angle derivative preserves order; output has same structure as input
        output_degrees =
            Vector{ResonantActionAnglePolynomialDegree}(undef, length(a.degrees))
        for (out_idx, degree) in enumerate(a.degrees)
            output_degrees[out_idx] = differentiate(degree, b)
        end
    end

    return ResonantActionAnglePolynomial(output_degrees)
end

"""
    Poisson Bracket
"""
function poisson_bracket(u::AbstractHomogenousPolynomial, v::AbstractHomogenousPolynomial)
    out_order = u.order+v.order-2
    out_polynomial = ZeroHomogenousPolynomial(out_order)
    for i = 1:2:5
        q = i
        p = i+1
        out_polynomial +=
            differentiate(u, q) * differentiate(v, p) -
            differentiate(u, p) * differentiate(v, q)
    end
    return out_polynomial
end


""" 
    Homogenous change in variables
"""
function changevars(
    a::AbstractHomogenousPolynomial,
    linchange_vars::Vector{<:AbstractHomogenousPolynomial},
)
    out_order = a.order+linchange_vars[1].order-1
    out_polynomial = ZeroHomogenousPolynomial(out_order)
    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)
    for i in eachindex(a_indices)
        a_multiindex = get_multiindex6(a.order, a_indices[i])
        out_polynomial += (
            a_coeffs[i] * prod([
                linchange_vars[1] ^ a_multiindex[1],
                linchange_vars[2] ^ a_multiindex[2],
                linchange_vars[3] ^ a_multiindex[3],
                linchange_vars[4] ^ a_multiindex[4],
                linchange_vars[5] ^ a_multiindex[5],
                linchange_vars[6] ^ a_multiindex[6],
            ])
        )
    end
    return out_polynomial
end

function changevars(
    a::MixedDegreePolynomial,
    linchange_vars::Vector{<:AbstractHomogenousPolynomial},
)
    output_degrees = AbstractHomogenousPolynomial[]
    for degree in a.degrees
        if nnz(degree.terms) != 0
            push!(output_degrees, changevars(degree, linchange_vars))
        end
    end
    return CombinePolynomials(output_degrees)
end

"""
    Substitution by number
"""
function evaluatepolynomial(a::RealHomogenousPolynomial, subs_vars::SVector{6})
    out_value = 0.0
    a_indices = rowvals(a.terms)
    a_coeffs = nonzeros(a.terms)
    for i in eachindex(a_indices)
        a_multiindex = get_multiindex6(a.order, a_indices[i])
        out_value+=a_coeffs[i]*prod(subs_vars .^ a_multiindex)
    end
    return out_value
end

function evaluatepolynomial(a::MixedDegreePolynomial, subs_vars::SVector{6})
    out_value = 0
    for degree in a.degrees
        if nnz(degree.terms) != 0
            out_value += evaluatepolynomial(degree, subs_vars)
        end
    end
    return out_value
end

function evaluatepolynomial(a::MixedDegreePolynomial, subs_vec::Vector{<:SVector})
    out_vec = zeros(eltype(subs_vec[1]), length(subs_vec))
    for degree in a.degrees
        if nnz(degree.terms) != 0
            degree_indices = rowvals(degree.terms)
            degree_coeffs = nonzeros(degree.terms)
            degree_order = degree.order
            for i in eachindex(degree_indices)
                degree_multiindex = get_multiindex6(degree_order, degree_indices[i])
                out_vec += [
                    degree_coeffs[i] * prod(subs_vars .^ degree_multiindex) for
                    subs_vars in subs_vec
                ]
            end
        end
    end
    return out_vec
end

function evaluatepolynomial(a::ResonantActionAnglePolynomialDegree, subs_vars::SVector{6})
    # Unpack subs vars 
    # get birkhoff coords from resonant
    I1 = subs_vars[1]
    I2 = subs_vars[3]
    I3 = subs_vars[5] - subs_vars[3]
    theta2 = subs_vars[4]
    action_vars = [I1, I2, I3]

    out_value = 0.0
    a_integrable_indices = rowvals(a.integrable_terms)
    a_integrable_coeffs = nonzeros(a.integrable_terms)
    a_resonant_indices = rowvals(a.resonant_terms)
    a_resonant_coeffs = nonzeros(a.resonant_terms)

    for i in eachindex(a_integrable_indices)
        a_multiindex = get_multiindex3(a.order, a_integrable_indices[i])
        out_value+=a_integrable_coeffs[i]*prod(action_vars .^ a_multiindex)
    end
    for i in eachindex(a_resonant_indices)
        a_resonant_coeff = a_resonant_coeffs[i]
        a_multiindex = get_multiindex3(a.order, a_resonant_indices[i])
        action_monomial = prod(action_vars .^ a_multiindex)
        if a_resonant_coeff[1] isa AbstractVector
            # SVector(num_coeffs::Vector, sin_freqs::Vector, cos_freqs::Vector)
            for j in eachindex(a_resonant_coeff[1])
                num = a_resonant_coeff[1][j]
                sin_freq = a_resonant_coeff[2][j]
                cos_freq = a_resonant_coeff[3][j]
                if sin_freq != 0
                    out_value += num * action_monomial * sin(sin_freq * theta2)
                elseif cos_freq != 0
                    out_value += num * action_monomial * cos(cos_freq * theta2)
                end
            end
        else
            if a_resonant_coeff[2] != 0
                out_value += a_resonant_coeff[1] * action_monomial * sin(a_resonant_coeff[2] * theta2)
            elseif a_resonant_coeff[3] != 0
                out_value += a_resonant_coeff[1] * action_monomial * cos(a_resonant_coeff[3] * theta2)
            end
        end
    end
    return out_value
end

function evaluatepolynomial(a::ResonantActionAnglePolynomial, subs_vars::SVector{6})
    out_value = 0
    for degree in a.degrees
        if nnz(degree.integrable_terms) != 0 || nnz(degree.resonant_terms) != 0
            out_value += evaluatepolynomial(degree, subs_vars)
        end
    end
    return out_value
end

function evaluatepolynomial(a::ResonantActionAnglePolynomial, subs_vec::Vector{<:SVector})
    out_vec = zeros(eltype(subs_vec[1]), length(subs_vec))
    for degree in a.degrees
        if nnz(degree.integrable_terms) != 0 || nnz(degree.resonant_terms) != 0
            degree_integrable_indices = rowvals(degree.integrable_terms)
            degree_integrable_coeffs = nonzeros(degree.integrable_terms)
            degree_resonant_indices = rowvals(degree.resonant_terms)
            degree_resonant_coeffs = nonzeros(degree.resonant_terms)
            for i in eachindex(degree_integrable_indices)
                degree_multiindex =
                    get_multiindex3(degree.order, degree_integrable_indices[i])
                out_vec += [
                    degree_integrable_coeffs[i]*prod(
                        [subs_vars[1], subs_vars[3], subs_vars[5]-subs_vars[3]] .^
                        degree_multiindex,
                    ) for subs_vars in subs_vec
                ]
            end
            for i in eachindex(degree_resonant_indices)
                degree_resonant_coeff = degree_resonant_coeffs[i]
                degree_multiindex =
                    get_multiindex3(degree.order, degree_resonant_indices[i])
                if degree_resonant_coeff[1] isa AbstractVector
                    for j in eachindex(degree_resonant_coeff[1])
                        num = degree_resonant_coeff[1][j]
                        sin_freq = degree_resonant_coeff[2][j]
                        cos_freq = degree_resonant_coeff[3][j]
                        if sin_freq != 0
                            out_vec += [
                                num * prod(
                                    [subs_vars[1], subs_vars[3], subs_vars[5]-subs_vars[3]] .^
                                    degree_multiindex,
                                ) * sin(sin_freq * subs_vars[4]) for subs_vars in subs_vec
                            ]
                        elseif cos_freq != 0
                            out_vec += [
                                num * prod(
                                    [subs_vars[1], subs_vars[3], subs_vars[5]-subs_vars[3]] .^
                                    degree_multiindex,
                                ) * cos(cos_freq * subs_vars[4]) for subs_vars in subs_vec
                            ]
                        end
                    end
                else
                    if degree_resonant_coeff[2] != 0
                        out_vec += [
                            degree_resonant_coeff[1]*prod(
                                [subs_vars[1], subs_vars[3], subs_vars[5]-subs_vars[3]] .^
                                degree_multiindex,
                            )*sin(degree_resonant_coeff[2]*subs_vars[4]) for
                            subs_vars in subs_vec
                        ]
                    elseif degree_resonant_coeff[3] != 0
                        out_vec += [
                            degree_resonant_coeff[1]*prod(
                                [subs_vars[1], subs_vars[3], subs_vars[5]-subs_vars[3]] .^
                                degree_multiindex,
                            )*cos(degree_resonant_coeff[3]*subs_vars[4]) for
                            subs_vars in subs_vec
                        ]
                    end
                end
            end
        end
    end
    return out_vec
end


function evaluatepolynomial(a::Vector{<:RealHomogenousPolynomial}, subs_vars::SVector{6})
    out_vec = Vector{eltype(subs_vars)}(undef, length(a))
    Threads.@threads for i = 1:length(a)
        out_vec[i] = evaluatepolynomial(a[i], subs_vars)
    end
    return out_vec
end

function evaluatepolynomial(a::Vector{<:MixedDegreePolynomial}, subs_vars::SVector{6})
    out_vec = Vector{eltype(subs_vars)}(undef, length(a))
    Threads.@threads for i = 1:length(a)
        out_vec[i] = evaluatepolynomial(a[i], subs_vars)
    end
    return out_vec
end

function evaluatepolynomial(a::Matrix{<:RealHomogenousPolynomial}, subs_vars::SVector{6})
    out_mat = Matrix{eltype(subs_vars)}(undef, size(a)[1], size(a)[2])
    for i = 1:size(a)[1]
        Threads.@threads for ii = 1:size(a)[2]
            out_mat[i, ii] = evaluatepolynomial(a[i, ii], subs_vars)
        end
    end
    return out_mat
end

function evaluatepolynomial(a::Matrix{<:MixedDegreePolynomial}, subs_vars::SVector{6})
    out_mat = Matrix{eltype(subs_vars)}(undef, size(a)[1], size(a)[2])
    for i = 1:size(a)[1]
        Threads.@threads for ii = 1:size(a)[2]
            out_mat[i, ii] = evaluatepolynomial(a[i, ii], subs_vars)
        end
    end
    return out_mat
end
