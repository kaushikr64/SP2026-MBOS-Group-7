function Base.show(io::IO, polynomial::AbstractHomogenousPolynomial)
    if nnz(polynomial.terms) == 1
        println(
            "Homogenous polynomial of order ",
            polynomial.order,
            " containing ",
            nnz(polynomial.terms),
            " term.",
        )
    else
        println(
            "Homogenous polynomial of order ",
            polynomial.order,
            " containing ",
            nnz(polynomial.terms),
            " terms.",
        )
    end
    polynomial_nzterms = findnz(polynomial.terms)
    polynomial_indices = polynomial_nzterms[1]
    polynomial_coeffs = polynomial_nzterms[2]
    for i = 1:lastindex(polynomial_indices)
        println(
            "Coefficient = ",
            polynomial_coeffs[i],
            ", Exponent = ",
            get_multiindex6(polynomial.order, polynomial_indices[i]),
        )
    end
end

function Base.show(io::IO, polynomial::MixedDegreePolynomial)
    if length(polynomial.degrees) == 1
        println(
            "Mixed degree polynomial containing ",
            length(polynomial.degrees),
            " homogenous polynomial.",
        )
    else
        println(
            "Mixed degree polynomial containing ",
            length(polynomial.degrees),
            " homogenous polynomials.",
        )
    end
    for degree in polynomial.degrees
        if nnz(degree.terms) != 0
            if nnz(degree.terms) == 1
                println("Order ", degree.order, ": ", nnz(degree.terms), " terms")
            else
                println("Order ", degree.order, ": ", nnz(degree.terms), " terms")
            end
        end
    end
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    polynomial_vector::AbstractVector{MixedDegreePolynomial},
)
    println("Vector with ", length(polynomial_vector), " polynomial elements:")
    for i in eachindex(polynomial_vector)
        polynomial = polynomial_vector[i]
        println("\nVector element ", i)
        show(polynomial)
    end
end


function Base.show(
    io::IO,
    ::MIME"text/plain",
    polynomial_matrix::AbstractMatrix{MixedDegreePolynomial},
)
    rows, cols = size(polynomial_matrix)
    println("Matrix with ", rows, "×", cols, " polynomial elements:")
    for i = 1:rows
        for j = 1:cols
            println(io, "\nMatrix element [", i, ",", j, "]")
            show(polynomial_matrix[i, j])
        end
    end
end


function Base.show(io::IO, polynomial::BirkhoffActionAnglePolynomial)
    println(
        "Birkhoff action-angle polynomial of degree ",
        length(polynomial.degrees) - 1,
    )
    for degree in polynomial.degrees
        if nnz(degree.terms) != 0
            println("Order ", degree.order, ": ", nnz(degree.terms), " terms")
            degree_nzterms = findnz(degree.terms)
            degree_indices = degree_nzterms[1]
            degree_coeffs = degree_nzterms[2]
            for i = 1:lastindex(degree_indices)
                println(
                    "  Coefficient = ",
                    degree_coeffs[i],
                    ", Exponent = ",
                    get_multiindex3(degree.order, degree_indices[i]),
                )
            end
        end
    end
end


function Base.show(io::IO, polynomial::ResonantActionAnglePolynomialDegree)
    println(
        "Resonant action-angle polynomial terms of degree ",
        polynomial.order,
        " containing ",
        nnz(polynomial.integrable_terms),
        " integrable terms and ",
        nnz(polynomial.resonant_terms),
        " resonant terms.",
    )
    integrable_nzterms = findnz(polynomial.integrable_terms)
    integrable_indices = integrable_nzterms[1]
    integrable_coeffs = integrable_nzterms[2]
    println("Integrable terms: ")
    for i = 1:lastindex(integrable_indices)
        println(
            "Coefficient = ",
            integrable_coeffs[i],
            ", Exponent = ",
            get_multiindex3(polynomial.order, integrable_indices[i]),
        )
    end
    if nnz(polynomial.resonant_terms) != 0
        resonant_nzterms = findnz(polynomial.resonant_terms)
        resonant_indices = resonant_nzterms[1]
        resonant_coeffs = resonant_nzterms[2]
        println("Resonant terms: ")
        for i = 1:lastindex(resonant_indices)
            exponent = get_multiindex3(polynomial.order, resonant_indices[i])
            coeff = resonant_coeffs[i]
            if coeff[1] isa AbstractVector
                # SVector(num_coeffs::Vector, trig1::Vector, trig2::Vector)
                for j = 1:length(coeff[1])
                    println(
                        "  Coefficient = ", coeff[1][j],
                        ", sin_freq = ", coeff[2][j],
                        ", cos_freq = ", coeff[3][j],
                        ", Exponent = ", exponent,
                    )
                end
            else
                println(
                    "  Coefficient = ", coeff[1],
                    ", sin_freq = ", coeff[2],
                    ", cos_freq = ", coeff[3],
                    ", Exponent = ", exponent,
                )
            end
        end
    end
end

function Base.show(io::IO, polynomial::ResonantActionAnglePolynomial)
    println("Resonant action-angle polynomial of degree ", length(polynomial.degrees)-1)
    for degree in polynomial.degrees
        if nnz(degree.integrable_terms) != 0 && nnz(degree.resonant_terms) != 0
            println(
                "Order ",
                degree.order,
                ": ",
                nnz(degree.integrable_terms),
                " integrable terms, and ",
                nnz(degree.resonant_terms),
                " resonant terms",
            )
        elseif nnz(degree.integrable_terms) != 0
            println(
                "Order ",
                degree.order,
                ": ",
                nnz(degree.integrable_terms),
                " integrable terms",
            )
        elseif nnz(degree.resonant_terms) != 0
            println(
                "Order ",
                degree.order,
                ": ",
                nnz(degree.resonant_terms),
                " resonant terms",
            )
        end
    end
end




function unit_variableset(length::Integer = 6)
    noexponent = zeros(Integer, length)
    variableset = Vector{RealHomogenousPolynomial}(undef, length)
    for i = 1:length
        exponent = copy(noexponent)
        exponent[i] = 1
        variableset[i] = RealHomogenousPolynomial(1, exponent);
    end
    return variableset
end

function terms_ofdegree(polynomial::MixedDegreePolynomial, degree::Integer)
    return polynomial.degrees[degree+1]
end

function terms_ofdegree(polynomial::ResonantActionAnglePolynomial, degree::Integer)
    return polynomial.degrees[degree+1]
end

function get_maxorder(polynomial::MixedDegreePolynomial)
    return findlast(polyterms -> nnz(polyterms.terms) != 0, polynomial.degrees) - 1
end
