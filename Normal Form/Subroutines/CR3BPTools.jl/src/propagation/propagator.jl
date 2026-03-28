"""
    Wrapper for propagation similar to what's there in MATLAB
"""
function propagate(f::Function,x0,tspan::Tuple,p;method = DP8(),reltol = 1e-12, abstol = 1e-12,kwargs...)
    prob = ODEProblem(f,x0,tspan,p)
    return solve(prob,method,reltol = reltol,abstol = abstol;kwargs...)
end

function propagate(f::Function,x0,tspan::Union{Vector,LinRange},p;method = DP8(),reltol = 1e-12, abstol = 1e-12,kwargs...)
    prob = ODEProblem(f,x0,[tspan[1],tspan[end]],p)
    sol = solve(prob,method,reltol = reltol,abstol = abstol;kwargs...)
    return sol(tspan).u
end

function propagate(f::Function,x0,tf,p;method = DP8(),reltol = 1e-12, abstol = 1e-12,kwargs...)
    prob = ODEProblem(f,x0,(0,tf),p)
    sol = solve(prob,method,reltol = reltol,abstol = abstol,kwargs...)
    return sol(tf)
end

function init_integrator(f::Function,x0,tspan::Tuple,p;method = DP8(),reltol = 1e-12, abstol = 1e-12,callback = nothing, kwargs...)
    prob = ODEProblem(f,x0,tspan,p)
    return init(prob,method,reltol = reltol,abstol = abstol;kwargs...)
end

function init_integrator(f::Function,x0,tf,p;method = DP8(),reltol = 1e-12, abstol = 1e-12,callback = nothing,kwargs...)
    prob = ODEProblem(f,x0,(0,tf),p)
    return init(prob,method,reltol = reltol,abstol = abstol,kwargs...)
end