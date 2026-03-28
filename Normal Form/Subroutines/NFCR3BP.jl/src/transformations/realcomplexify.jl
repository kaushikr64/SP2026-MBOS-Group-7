# The complex change changes the order of variables
function complexifying_change()
    return 1/sqrt(2) * [
        sqrt(2) 0 0 0 0 0
        0 0 1 1im 0 0
        0 0 0 0 1 1im
        0 sqrt(2) 0 0 0 0
        0 0 1im 1 0 0
        0 0 0 0 1im 1
    ]
end

# The realifying change does not change the order of variables
function realifying_change()
    return 1/sqrt(2) * [
        sqrt(2) 0 0 0 0 0
        0 sqrt(2) 0 0 0 0
        0 0 1 -1im 0 0
        0 0 -1im 1 0 0
        0 0 0 0 1 -1im
        0 0 0 0 -1im 1
    ]
end

function realify(a::ComplexHomogenousPolynomial)
    realifying_transformation = realifying_change()*unit_variableset()
    return RealHomogenousPolynomial(
        a.order,
        real.(changevars(a, realifying_transformation).terms),
    )
end

function realify(a::MixedDegreePolynomial)
    realifying_transformation = realifying_change()*unit_variableset()
    c_degrees = Vector{RealHomogenousPolynomial}(undef, length(a.degrees))
    for i in eachindex(a.degrees)
        order = i-1
        degree = a.degrees[i]
        c_degrees[i] = RealHomogenousPolynomial(
            order,
            real.(changevars(degree, realifying_transformation).terms),
        )
    end
    return CombinePolynomials(c_degrees)
end
