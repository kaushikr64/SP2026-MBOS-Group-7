"""
    Evaluate Hamilton's equations in the order of coefficients qdot1,pdot1,qdot2,pdot2 etc
"""
function get_hamiltonian_flow(H::RealHomogenousPolynomial)
    dxdt = Vector{RealHomogenousPolynomial}(undef, 6)
    for i = 1:2:5
        dxdt[i] = differentiate(H, i+1)
        dxdt[i+1] = -1*differentiate(H, i+1)
    end
    return dxdt
end

function get_hamiltonian_flow(H::MixedDegreePolynomial)
    dxdt = Vector{MixedDegreePolynomial}(undef, 6)

    #qdot
    dxdt[1] = differentiate(H, 2)
    dxdt[3] = differentiate(H, 4)
    dxdt[5] = differentiate(H, 6)

    dxdt[2] = -1*differentiate(H, 1)
    dxdt[4] = -1*differentiate(H, 3)
    dxdt[6] = -1*differentiate(H, 5)

    return dxdt
end

function get_hamiltonian_jacobian(Hvec::Vector{MixedDegreePolynomial})
    jacbn = Matrix{MixedDegreePolynomial}(undef, 6, 6)
    for i = 1:6
        for ii = 1:6
            jacbn[i, ii] = differentiate(Hvec[i], ii)
        end
    end
    return jacbn
end

function get_hamiltonian_flow(H::ResonantActionAnglePolynomial)
    zero_resonantAApolynomial = ResonantActionAnglePolynomial([
        ResonantActionAnglePolynomialDegree(0, spzeros(Real, 1), spzeros(SVector{3}, 1)),
    ])

    # Action angle rates of change (angle flow)
    dHdI1 = differentiate(H, 1)
    dHdI2 = differentiate(H, 2)
    dHdI3 = differentiate(H, 3)
    dHdθ̂2 = differentiate(H, 4)

    dÎ1dt = zero_resonantAApolynomial
    dÎ2dt = -1*dHdθ̂2
    dÎ3dt = zero_resonantAApolynomial

    dθ̂1dt = dHdI1
    dθ̂2dt = dHdI2-dHdI3
    dθ̂3dt = dHdI3

    return [dÎ1dt, dθ̂1dt, dÎ2dt, dθ̂2dt, dÎ3dt, dθ̂3dt]
end
