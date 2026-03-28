"""
utilities related to the lie series transformation
"""
function birkhoff_term(k::AbstractVector{<:Integer})
    #Birkhoff
    if (k[1] == k[2]) && (k[3] == k[4]) && (k[5] == k[6])
        return true
    end
    return false
end

function resonant_term(k::AbstractVector{<:Integer})
    #Resonant
    if (k[1] == k[2]) && (k[3] - k[4] + k[5] - k[6] == 0)
        return true
    end
    return false
end


function obtain_generating_function(
    H::MixedDegreePolynomial,
    order::Integer,
    sys::NamedTuple,
)
    degreeindex = order+1
    eta = sys.eigenvalues
    remove_term = sys.removal_criteria
    h_terms = findnz(H.degrees[degreeindex].terms)
    h_indices = h_terms[1]
    h_coeffs = h_terms[2]
    g_terms = spzeros(Complex, get_num_monomials(6, order))
    indices_to_remove = []
    for i in eachindex(h_indices)
        listindex = h_indices[i]
        multiindex = get_multiindex6(order, listindex)
        if remove_term(multiindex)
            kq = multiindex[1:2:5]
            kp = multiindex[2:2:6]
            h_coeff = h_coeffs[i]
            g_coeff = -h_coeff/dot(kp-kq, eta)
            g_terms[listindex] = g_coeff
            push!(indices_to_remove, listindex)
        end
    end
    dropzeros!(H.degrees[degreeindex].terms)
    return ComplexHomogenousPolynomial(order, g_terms), indices_to_remove
end

function lie_series_transformation!(H::MixedDegreePolynomial, sys::NamedTuple)
    maxorder = sys.maxorder_NF
    G_list = AbstractHomogenousPolynomial[]
    for order = 3:maxorder
        println("Normalzing order ", order)
        generating_function, indices_to_remove = obtain_generating_function(H, order, sys)
        order_increment = order-2 #Increment in order resulting from poisson brackets
        for i = maxorder:-1:order
            degreeindex = i + 1
            stepstaken = maxorder - i #number of steps already taken for this order
            substepnum = div(stepstaken, order_increment) + 1 #number of substeps needed
            poiss_brac = terms_ofdegree(H, i-order_increment) #initialize poisson bracket recurrence
            for ii = 1:substepnum
                poiss_brac = poisson_bracket(poiss_brac, generating_function) # Do the poisson brackets
                H.degrees[degreeindex+(ii-1)*order_increment] +=
                    ((1/factorial(big(ii)))*poiss_brac)
                # Display counter
                # println("Step ", stepstaken+1, ".",ii)
            end
        end
        push!(G_list, generating_function)
        #Remove all terms predicted to be removed by the term removal criterion
        terms_ofdegree(H, order).terms[indices_to_remove].=0
        dropzeros!(terms_ofdegree(H, order).terms)
    end
    return H, CombinePolynomials(G_list)
end
