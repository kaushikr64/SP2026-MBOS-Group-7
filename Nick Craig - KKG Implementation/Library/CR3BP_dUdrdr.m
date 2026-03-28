function dUdrdr = CR3BP_dUdrdr(X,mu)
% This function returns the second derivative of the pseudo-potential
% function with respect to the position vector evaluated at a given
% position.
% 
% Inputs: 
%   r [x,y,z]'
%       x: x-coordinate in rotating frame
%       y: y-coordinate in rotating frame
%       z: z-coordinate in rotating frame
%   mu: 
%
% Outputs:
%   dUdrdr: pseudo-potential double derivative matrix
% 
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

if ~exist('mu','var')
    mu = 1.215059e-2; % Earth-Moon MU
end

x = X(1,1);
y = X(2,1);
z = X(3,1);

%{
The commented section below analytically computes the double derivative of
the pseudo-potential function. The results were hard coded into the
function to resduce computational complexity.

dUdrdr = hessian(CR3BP_U,[sym('x'),sym('y'),sym('z')]);
dUdrdr = double(dUdrdr(x,y,z,mu));
%}

dUdxdx = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2)...
    - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2)...
    + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2))...
    - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
dUdxdy = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2))...
    - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
dUdxdz = (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2))...
    - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));

dUdydy = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2)...
    - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2)...
    - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2)...
    + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
dUdydz = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)...
    - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);

dUdzdz = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2)...
    - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2)...
    - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2)...
    + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);

dUdrdr = [dUdxdx,dUdxdy,dUdxdz;...
    dUdxdy,dUdydy,dUdydz;...
    dUdxdz,dUdydz,dUdzdz]; 
end

