"""
    Action Angle transformations
"""

function normalform2birkhoff(complex_normal_form::MixedDegreePolynomial)
    maxorder_nf = get_maxorder(complex_normal_form)
    maxorder_aa = Int(floor(maxorder_nf / 2))
    birkhoff_degrees = Vector{BirkhoffActionAnglePolynomialDegree}(undef, maxorder_aa + 1)

    birkhoff_degrees[1] = BirkhoffActionAnglePolynomialDegree(0, spzeros(Float64, 1))

    for actionangle_degree = 1:maxorder_aa
        normalform_degree = 2 * actionangle_degree
        normalform_ofdegree = terms_ofdegree(complex_normal_form, normalform_degree)
        normalform_indices = rowvals(normalform_ofdegree.terms)
        normalform_coeffs = nonzeros(normalform_ofdegree.terms)

        num_monomials = get_num_monomials(3, actionangle_degree)
        terms = spzeros(Float64, num_monomials)

        for i in eachindex(normalform_indices)
            normalform_exponent = get_multiindex6(normalform_degree, normalform_indices[i])
            (kq1, _, kq2, kp2, kq3, kp3) = normalform_exponent
            actionangle_exponent = SVector(kq1, kq2, kq3)
            actionangle_index = get_listindex3(actionangle_degree, actionangle_exponent)
            terms[actionangle_index] = real(normalform_coeffs[i] / (1im^kp2 * 1im^kp3))
        end

        birkhoff_degrees[actionangle_degree + 1] = BirkhoffActionAnglePolynomialDegree(
            actionangle_degree, droptol!(terms, POLYTOL)
        )
    end

    return BirkhoffActionAnglePolynomial(birkhoff_degrees)
end

function normalform2birkhoff(normalform_state::AbstractVector)
    birkhoff_state = Vector{eltype(normalform_state)}(undef, 6)
    # Compute actions
    birkhoff_state[1] = normalform_state[1]
    birkhoff_state[3] = 0.5*(normalform_state[3]^2+normalform_state[4]^2)
    birkhoff_state[5] = 0.5*(normalform_state[5]^2+normalform_state[6]^2)

    # Compute angles
    birkhoff_state[2] = normalform_state[2]
    birkhoff_state[4] = -atan(normalform_state[4], normalform_state[3])
    birkhoff_state[6] = -atan(normalform_state[6], normalform_state[5])
    return SVector{6}(birkhoff_state)
end

function birkhoff2normalform(birkhoff_state::AbstractVector)
    normalform_state = Vector{eltype(birkhoff_state)}(undef, 6)

    normalform_state[1] = birkhoff_state[1]
    normalform_state[2] = birkhoff_state[2]
    normalform_state[3] = sqrt(2*birkhoff_state[3])*cos(birkhoff_state[4])
    normalform_state[4] = -sqrt(2*birkhoff_state[3])*sin(birkhoff_state[4])
    normalform_state[5] = sqrt(2*birkhoff_state[5])*cos(birkhoff_state[6])
    normalform_state[6] = -sqrt(2*birkhoff_state[5])*sin(birkhoff_state[6])
    return SVector{6}(normalform_state)
end

"""
    resonant
"""
function normalform2resonant(complex_normal_form::MixedDegreePolynomial)
    # Get the max order of resonant normal form
    maxorder_resnf = Int(floor(get_maxorder(complex_normal_form)/2))
    resonant_AA_degrees =
        Vector{ResonantActionAnglePolynomialDegree}(undef, maxorder_resnf+1)

    resonant_AA_degrees[1] =
        ResonantActionAnglePolynomialDegree(0, spzeros(Real, 1), spzeros(SVector{3}, 1))

    for actionangle_degree = 1:maxorder_resnf
        normalform_degree = 2*actionangle_degree
        normalform_ofdegree = terms_ofdegree(complex_normal_form, normalform_degree)
        normalform_indices = rowvals(normalform_ofdegree.terms)
        normalform_coeffs = nonzeros(normalform_ofdegree.terms)

        num_monomials = get_num_monomials(3, actionangle_degree)
        integrable_terms = spzeros(Real, num_monomials)
        resonant_terms = spzeros(SVector{3}, num_monomials)

        for i in eachindex(normalform_indices)
            normalform_exponent = get_multiindex6(normalform_degree, normalform_indices[i])
            (kq1, kp1, kq2, kp2, kq3, kp3) = normalform_exponent

            actionangle_exponent = SVector(kq1, Int((kq2+kp2)/2), Int((kq3+kp3)/2))
            actionangle_index = get_listindex3(actionangle_degree, actionangle_exponent)
            actionangle_num_coeff = real(normalform_coeffs[i]/(1im^kp2*1im^kp3))

            if kq2 == kp2 && kq3 == kp3
                integrable_terms[actionangle_index] = actionangle_num_coeff

            elseif kq2 > kp2 && kq2-kp2 == kp3-kq3
                if resonant_terms[actionangle_index][1] == 0
                    actionangle_coeff = SVector(2*actionangle_num_coeff, 0, kq2-kp2)
                elseif resonant_terms[actionangle_index][1] != 0
                    existing_actionangle_coeff = resonant_terms[actionangle_index]
                    actionangle_coeff = SVector(
                        vcat(existing_actionangle_coeff[1], 2*actionangle_num_coeff),
                        vcat(existing_actionangle_coeff[2], 0),
                        vcat(existing_actionangle_coeff[3], kq2-kp2),
                    )
                end
                resonant_terms[actionangle_index] = actionangle_coeff
            end
        end
        resonant_AA_degree = ResonantActionAnglePolynomialDegree(
            actionangle_degree,
            integrable_terms,
            resonant_terms,
        )
        resonant_AA_degrees[actionangle_degree+1] = resonant_AA_degree
    end

    return ResonantActionAnglePolynomial(resonant_AA_degrees)
end

function normalform2resonant(normalform_state::AbstractVector)
    birkhoff_state = normalform2birkhoff(normalform_state)

    resonant_state = @SVector [
        birkhoff_state[1],
        birkhoff_state[2],
        birkhoff_state[3],
        mod(birkhoff_state[4]-birkhoff_state[6], 2π),
        birkhoff_state[5]+birkhoff_state[3],
        mod(birkhoff_state[6], 2π),
    ]

    return resonant_state
end

function resonant2normalform(resonant_state::AbstractVector)
    birkhoff_state = @SVector [
        resonant_state[1]
        resonant_state[2]
        resonant_state[3]
        resonant_state[4]+resonant_state[6]
        resonant_state[5]-resonant_state[3]
        resonant_state[6]
    ]

    return birkhoff2normalform(birkhoff_state)
end
