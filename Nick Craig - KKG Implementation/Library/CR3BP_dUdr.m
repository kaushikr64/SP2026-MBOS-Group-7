function dUdr = CR3BP_dUdr(r,mu)
% This function returns the derivative of the pseudo-potential function
% with respect to the position vector evaluated at the position vector.
%
% Inputs: 
%   r [x,y,z]'
%       x: x-coordinate in rotating frame
%       y: y-coordinate in rotating frame
%       z: z-coordinate in rotating frame
%   mu: 
%
% Outputs:
%   dUdr: pseudo-potential derivative vector
% 
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

if ~exist('mu','var')
    mu = 1.215059e-2; % Earth-Moon MU
end

x = r(1,1);
y = r(2,1);
z = r(3,1);

%{
The commented section below analytically computes the derivative of
the pseudo-potential function. The results were hard coded into the
function to resuce computational complexity.

dUdr = jacobian(CR3BP_U,[sym('x'),sym('y'),sym('z')]);
dUdr = double(dUdr(x,y,z,mu))';
%}

dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2))...
    - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));

dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2)...
    + (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);

dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2)...
    - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);

dUdr = [dUdx;dUdy;dUdz];
end

