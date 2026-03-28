function X_dot = CR3BP_STM_EOM(X,mu,U)
% Author: Nick Craig
% Date: 03/31/2025
% 
% Description:
% This function computes the numerical equations of motion for a body
% acting within the Circular Restriced 3 Body Problem. Here, we also
% include the state transition matrix, and if requested, the linearized
% control B matrix.
%
% Inputs: 
%   X       [42,1]  Includes vectorized STM
%           [60,1]  Includes vectorized STM and B  
%   mu      [1]     CR3BP non-dimensional system gravitational parameter
%   U       [3,1]   Object external acceleration
% 
% Ouputs: 
%   X_dot   [42,1]  Includes vectorized STM derivative
%           [60,1]  Includes vectorized B derivative
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

if ~exist('U','var')
    U = [0;0;0]; % Small mass external acceleration
end

nx = 6;
nu = 3;

r = X(1:3,1); % SC position
v = X(4:6,1); % SC velocity


rv_dot = CR3BP_EOM([r;v],mu,U);

if length(X) == 42
    STM = reshape(X(nx+1:nx+nx^2,1),[nx,nx]);
    dfdx = [zeros(3,3),eye(3);CR3BP_dUdrdr(r,mu),[0,2,0;-2,0,0;0,0,0]];
    
    STM_dot = dfdx*STM;

    X_dot = [rv_dot;reshape(STM_dot,[nx^2,1])];
end
if length(X) == 60
    STM = reshape(X(nx+1:nx+nx^2,1),[nx,nx]);
    dfdx = [zeros(3,3),eye(3);CR3BP_dUdrdr(r,mu),[0,2,0;-2,0,0;0,0,0]];
    STM_dot = dfdx*STM;
    
    B = [zeros(nu,nu);eye(nu)];
    Bdot = STM\B;

    X_dot = [rv_dot;reshape(STM_dot,[nx^2,1]);reshape(Bdot,[nx*nu,1])];
end
end