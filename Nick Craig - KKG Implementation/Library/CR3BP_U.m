function U = CR3BP_U()
% This function returns the pseudo-potential function for the CR3BP system
% in symbolic form. There are no inputs and the output is a symbolic
% funciton with 4 inputs
% 
% Inputs:
%   N/A
% 
% Outputs:
%   Us(x,y,z,mu): symbolic pseudo-potential function
%       x:  x-coordinate of SC in rotating CR3BP frame
%       y:  y-coordinate of SC in rotating CR3BP frame
%       z:  z-coordinate of SC in rotating CR3BP frame
%       mu: non-dim constant 
% 
% ------------------------------------------------------------------------

syms x y z mu
U(x,y,z,mu) = (1-mu)/(((x+mu)^2+y^2+z^2)^(1/2))+mu/(((x-1+mu)^2+y^2+z^2)^(1/2))+...
    1/2*(x^2+y^2);
end
