function JC = CR3BP_JC(X,mu)
% This function returns the Jacobian Constant of a spacecraft at a point in
% the CR3BP system.
%
% Inputs: 
%   rv [x,y,z,xdot,ydot,zdot]'
%       x: x-coordinate in rotating frame
%       y: y-coordinate in rotating frame
%       z: z-coordinate in rotating frame
%       xdot: x-velocity in rotating frame
%       ydot: y-velocity in rotating frame
%       zdot: z-velocity in rotating frame
%   mu: 
%
% Outputs:
%   JC: Jacobian Constant
% 
% ------------------------------------------------------------------------

% JC = 2*double(subs(CR3BP_U,[sym('x'),sym('y'),sym('z'),sym('mu')],...
%     [X(1),X(2),X(3),mu]))-X(4:6)'*X(4:6);

x = X(1,1);
y = X(2,1);
z = X(3,1);
xdot = X(4,1);
ydot = X(5,1);
zdot = X(6,1);

d = sqrt((x+mu)^2+y^2+z^2);
r = sqrt((x-1+mu)^2+y^2+z^2);

vsqrd = xdot^2+ydot^2+zdot^2;

JC = x^2+y^2+2*(1-mu)/d + 2*mu/r -vsqrd;
end

