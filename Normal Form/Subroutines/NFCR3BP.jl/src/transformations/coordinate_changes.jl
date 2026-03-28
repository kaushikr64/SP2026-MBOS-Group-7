""" All changes in coordinates """



# Synodic -> Normalized direction
"""
    Synodic to recentered and vice versa
"""
function synodic2recentered(synodic_state::SVector{6}, transformations::NamedTuple)
    b = transformations.recentering
    T = transformations.rescaling
    V = transformations.qqdot2qp

    return SVector(inv(T)*(inv(V)*synodic_state-b)...)
end

function synodic2recentered(
    synodic_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    b = transformations.recentering
    T = transformations.rescaling
    V = transformations.qqdot2qp

    recentered_statehist =
        Vector{eltype(synodic_statehist)}(undef, length(synodic_statehist))
    for i in eachindex(synodic_statehist)
        recentered_statehist[i] = inv(T)*(inv(V)*synodic_statehist[i]-b)
    end
    return recentered_statehist
end

function recentered2synodic(recentered_state::SVector{6}, transformations::NamedTuple)
    b = transformations.recentering
    T = transformations.rescaling
    V = transformations.qqdot2qp

    return SVector{6}(V*(T*recentered_state+b))
end

function recentered2synodic(
    recentered_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    b = transformations.recentering
    T = transformations.rescaling
    V = transformations.qqdot2qp

    synodic_statehist =
        Vector{eltype(recentered_statehist)}(undef, length(recentered_statehist))
    for i in eachindex(synodic_statehist)
        synodic_statehist[i] = V*(T*recentered_statehist[i]+b)
    end
    return synodic_statehist
end

"""
    recentered to diagonal and vice verswa
"""
function recentered2diagonal(recentered_state::SVector{6}, transformations::NamedTuple)
    C_tf = transformations.diagonalizing

    return SVector(inv(C_tf)*recentered_state...)
end

function recentered2diagonal(
    recentered_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    C_tf = transformations.diagonalizing

    diagonal_statehist =
        Vector{eltype(recentered_statehist)}(undef, length(recentered_statehist))
    for i in eachindex(recentered_statehist)
        diagonal_statehist[i] = inv(C_tf)*recentered_statehist[i]
    end
    return diagonal_statehist
end

function diagonal2recentered(diagonal_state::SVector{6}, transformations::NamedTuple)
    C_tf = transformations.diagonalizing

    return SVector{6}(C_tf*diagonal_state)
end

function diagonal2recentered(
    diagonal_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    C_tf = transformations.diagonalizing

    recentered_statehist =
        Vector{eltype(diagonal_statehist)}(undef, length(diagonal_statehist))
    for i in eachindex(diagonal_statehist)
        recentered_statehist[i] = C_tf*diagonal_statehist[i]
    end
    return recentered_statehist
end

"""
    diagonal to normal form
"""
function diagonal2normalform(diagonal_state::SVector{6}, transformations::NamedTuple)
    qp_state = SVector(
        diagonal_state[1],
        diagonal_state[4],
        diagonal_state[2],
        diagonal_state[5],
        diagonal_state[3],
        diagonal_state[6],
    )

    backward_transformations = transformations.backward

    nf_hist = Vector{eltype(diagonal_state)}(undef, 6)

    Threads.@threads for i = 1:6
        nf_hist[i] = evaluatepolynomial(backward_transformations[i], qp_state)
    end

    normalform_state = SVector{6}(nf_hist)

    return normalform_state
end

function diagonal2normalform(
    diagonal_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    qp_statehist = [
        SVector(
            diagonal_vec[1],
            diagonal_vec[4],
            diagonal_vec[2],
            diagonal_vec[5],
            diagonal_vec[3],
            diagonal_vec[6],
        ) for diagonal_vec in diagonal_statehist
    ]

    backward_transformations = transformations.backward

    nf_hist = Vector{Vector}(undef, 6)

    Threads.@threads for i = 1:6
        nf_hist[i] = evaluatepolynomial(backward_transformations[i], qp_statehist)
    end

    normalform_statehist = [
        SVector(x, y, z, px, py, pz) for (x, y, z, px, py, pz) in
        zip(nf_hist[1], nf_hist[2], nf_hist[3], nf_hist[4], nf_hist[5], nf_hist[6])
    ]

    return normalform_statehist
end

function diagonal2normalform_numerical(
    diagonal_state::SVector{6},
    transformations::NamedTuple,
)
    G_flow = transformations.G_flow
    flowed_state = SVector(
        diagonal_state[1],
        diagonal_state[4],
        diagonal_state[2],
        diagonal_state[5],
        diagonal_state[3],
        diagonal_state[6],
    )

    for i = 3:MAXDEGREE
        if nnz(G_flow[1].degrees[i].terms) != 0
            flowed_state = propagate(
                eoms_polynomial,
                flowed_state,
                -1,
                [dgdt.degrees[i] for dgdt in G_flow],
            )
        end
    end
    return flowed_state
end

function diagonal2normalform_numerical(
    diagonal_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    qp_statehist = [
        SVector(
            diagonal_vec[1],
            diagonal_vec[4],
            diagonal_vec[2],
            diagonal_vec[5],
            diagonal_vec[3],
            diagonal_vec[6],
        ) for diagonal_vec in diagonal_statehist
    ]

    G_flow = transformations.G_flow

    normalform_statehist =
        Vector{eltype(diagonal_statehist)}(undef, length(diagonal_statehist))

    for i = 1:length(diagonal_statehist)
        flowed_state = qp_statehist[i]
        for i = 3:MAXDEGREE
            if nnz(G_flow[1].degrees[i].terms) != 0
                flowed_state = propagate(
                    eoms_polynomial,
                    flowed_state,
                    -1,
                    [dgdt.degrees[i] for dgdt in G_flow],
                )
            end
        end
        normalform_statehist[i] = flowed_state
    end

    return normalform_statehist
end


function normalform2diagonal(normalform_state::SVector{6}, transformations::NamedTuple)
    forward_transformations = transformations.forward

    qp_state = Vector{eltype(normalform_state)}(undef, 6)

    Threads.@threads for i = 1:6
        qp_state[i] = evaluatepolynomial(forward_transformations[i], normalform_state)
    end

    diagonal_statehist = SVector(
        qp_state[1],
        qp_state[3],
        qp_state[5],
        qp_state[2],
        qp_state[4],
        qp_state[6],
    )

    return diagonal_statehist
end

function normalform2diagonal(
    normalform_statehist::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    forward_transformations = transformations.forward

    qp_hist = Vector{Vector}(undef, 6)

    Threads.@threads for i = 1:6
        qp_hist[i] = evaluatepolynomial(forward_transformations[i], normalform_statehist)
    end

    diagonal_statehist = [
        SVector(x, y, z, px, py, pz) for (x, px, y, py, z, pz) in
        zip(qp_hist[1], qp_hist[2], qp_hist[3], qp_hist[4], qp_hist[5], qp_hist[6])
    ]

    return diagonal_statehist
end

"""
    normal form to action angle
"""
function normalform2actionangle(normalform_state::SVector{6}, transformations::NamedTuple)
    nf2aa = transformations.nf2aa
    return nf2aa(normalform_state)
end

function normalform2actionangle(
    normalform_state::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    nf2aa = transformations.nf2aa
    return nf2aa.(normalform_state)
end

function actionangle2normalform(actionangle_state::SVector{6}, transformations::NamedTuple)
    aa2nf = transformations.aa2nf
    return aa2nf(actionangle_state)
end

function actionangle2normalform(
    actionangle_state::AbstractVector{<:SVector},
    transformations::NamedTuple,
)
    aa2nf = transformations.aa2nf
    return aa2nf.(actionangle_state)
end

"""
    convert from pure action angles to working action angles
"""
function actionangle2pureactionangle(actionangle_state::SVector{6})
    return vcat(
        actionangle_state[1]*actionangle_state[2],
        -log(sqrt(normalform_state[1]/normalform_state[2])),
        actionangle_state[3:4]...,
    )
end

function actionangle2pureactionangle(actionangle_state::AbstractVector{<:SVector})
    return actionangle2pureactionangle.(actionangle_state)
end

function pureactionangle2actionangle(pureactionangle_state::SVector{6})
    return vcat(
        sqrt(pureactionangle_state[1])*exp(pureactionangle_state[2]),
        sqrt(pureactionangle_state[1])*exp(-pureactionangle_state[2]),
        pureactionangle_state[3:4]...,
    )
end

function pureactionangle2actionangle(pureactionangle_state::AbstractVector{<:SVector})
    return pureactionangle2actionangle.(pureactionangle_state)
end


"""synodic to normal form, action angle and vice versa. These are our most used transformations 
for stochastic and deterministic respectively"""
function synodic2normalform(synodic_stateorvec, transformations::NamedTuple)
    return diagonal2normalform(
        recentered2diagonal(
            synodic2recentered(synodic_stateorvec, transformations),
            transformations,
        ),
        transformations,
    )
end

function synodic2normalform_numerical(synodic_stateorvec, transformations::NamedTuple)
    return diagonal2normalform_numerical(
        recentered2diagonal(
            synodic2recentered(synodic_stateorvec, transformations),
            transformations,
        ),
        transformations,
    )
end

function normalform2synodic(normalform_stateorvec, transformations::NamedTuple)
    return recentered2synodic(
        diagonal2recentered(
            normalform2diagonal(normalform_stateorvec, transformations),
            transformations,
        ),
        transformations,
    )
end

function synodic2actionangle(synodic_stateorvec, transformations::NamedTuple)
    return normalform2actionangle(
        synodic2normalform(synodic_stateorvec, transformations),
        transformations,
    )
end

function synodic2actionangle_numerical(synodic_stateorvec, transformations::NamedTuple)
    return normalform2actionangle(
        synodic2normalform_numerical(synodic_stateorvec, transformations),
        transformations,
    )
end

function actionangle2synodic(actionangle_stateorvec, transformations::NamedTuple)
    return normalform2synodic(
        actionangle2normalform(actionangle_stateorvec, transformations),
        transformations,
    )
end
