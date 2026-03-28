"""
    equations of motion
"""
function eoms_CR3BP(rv,p::NamedTuple,t)    
    mu = p.mu
    x,y,z,xdot,ydot,zdot = rv

    #Vector from Mass 1
    r1 = norm([x+mu, y, z])

    #Vector from Mass 2
    r2 = norm([x-1+mu, y, z])


    #Calculating accelerations
    xddot = -(1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3 + 2*ydot + x
    yddot = -(1-mu)*y/r1^3 - mu*y/r2^3 - 2*xdot + y
    zddot = -(1-mu)*z/r1^3 - mu*z/r2^3
    

    return SVector(xdot,ydot,zdot,xddot,yddot,zddot)
end

function eoms_CR3BP(rv,p::Real,t)    
    mu = p
    x,y,z,xdot,ydot,zdot = rv

    #Vector from Mass 1
    r1 = norm([x+mu, y, z])

    #Vector from Mass 2
    r2 = norm([x-1+mu, y, z])


    #Calculating accelerations
    xddot = -(1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3 + 2*ydot + x
    yddot = -(1-mu)*y/r1^3 - mu*y/r2^3 - 2*xdot + y
    zddot = -(1-mu)*z/r1^3 - mu*z/r2^3
    

    return SVector(xdot,ydot,zdot,xddot,yddot,zddot)
end
