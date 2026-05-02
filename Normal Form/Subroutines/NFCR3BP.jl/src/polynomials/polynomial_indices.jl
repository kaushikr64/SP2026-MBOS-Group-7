#Psi Calculation table
function psi_calc(i::Integer, n::Integer)
    return binomial(n+i-1, i-1)
end

function psi_table(maxorder)
    table = Matrix{Integer}(undef, 5, maxorder+1)
    for i = 2:6
        for n = 0:maxorder
            table[i-1, n+1] = psi_calc(i, n)
        end
    end
    return table
end

const PSI = psi_table(64)

function get_num_monomials(numcoords, order)
    return (PSI[numcoords-1, order .+ 1])
end

#MultiIndex table
function successive_multiindex(exponents::AbstractVector{<:Integer})
    multiindex = collect(exponents)
    if exponents[1] != 0
        multiindex[1] -= 1
        multiindex[2] += 1
        return SVector(multiindex...)
    end
    if exponents[2] != 0
        multiindex[1] = multiindex[2]-1
        multiindex[2] = 0
        multiindex[3] += 1
        return SVector(multiindex...)
    end
    if exponents[3] != 0
        multiindex[1] = multiindex[3]-1
        multiindex[3] = 0
        multiindex[4] += 1
        return SVector(multiindex...)
    end
    if exponents[4] != 0
        multiindex[1] = multiindex[4]-1
        multiindex[4] = 0
        multiindex[5] += 1
        return SVector(multiindex...)
    end
    if exponents[5] != 0
        multiindex[1] = multiindex[5]-1
        multiindex[5] = 0
        multiindex[6] += 1
        return SVector(multiindex...)
    end
end

function generate_multiindex_table(maxorder::Integer, numcoords)
    table = Vector(undef, maxorder+1)
    for i = 0:maxorder
        num_monomials = get_num_monomials(numcoords, i)
        ordervec = Vector(undef, num_monomials)
        ordervec[1] = SVector(i, zeros(Int, numcoords-1)...)
        for ii = 2:num_monomials
            ordervec[ii] = successive_multiindex(ordervec[ii-1])
        end
        table[i+1] = ordervec
    end
    return table
end

const MAXDEGREE = 16
const MIDX6 = generate_multiindex_table(MAXDEGREE, 6)
const MIDX3 = generate_multiindex_table(MAXDEGREE, 3)


function get_multiindex6(order::Integer, listindex::Integer)
    return MIDX6[order+1][listindex]::SVector{6,Int}
end

function get_listindex6(order::Integer, multiindex)
    numcoords = 6
    listindex = 1
    n_i_prev = order
    for i = numcoords:-1:3
        # number of monomials of the indicies other than that considered
        k_i = multiindex[i]
        n_i = n_i_prev-k_i

        # Get the contribution for the current order
        listindex += sum(get_num_monomials(i-1, (n_i+1):n_i_prev))

        n_i_prev = n_i
    end
    # contribution from final term
    listindex += multiindex[2]
end

function get_multiindex3(order::Integer, listindex::Integer)
    return MIDX3[order+1][listindex]::SVector{3,Int}
end

function get_listindex3(order::Integer, multiindex)
    numcoords = 3
    listindex = 1
    n_i_prev = order
    for i = numcoords:-1:3
        # number of monomials of the indicies other than that considered
        k_i = multiindex[i]
        n_i = n_i_prev-k_i

        # Get the contribution for the current order
        listindex += sum(get_num_monomials(i-1, (n_i+1):n_i_prev))

        n_i_prev = n_i
    end
    # contribution from final term
    listindex += multiindex[2]
end
