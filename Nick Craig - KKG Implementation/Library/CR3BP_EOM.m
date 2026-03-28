function X_dot = CR3BP_EOM(X,mu,U)
% Author: Nick Craig
% Date: 03/31/2025
% 
% Description:
% This function initializes classic equations of motion for the state of a
% body acting within the Circular Restricted 3 Body Problem (CR3BP). 
%
% Inputs: 
%   X       [6,1]   Object state (r,v) in CR3BP rotating frame      
%   mu      [1]     CR3BP non-dimensional system gravitational parameter
%   U       [3,1]   Object external acceleration
% 
% Ouputs: 
%   Xdot    [6,1]   Object state time derivative in CR3BP rotating frame
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

if ~exist('U','var')
    U = [0;0;0]; % Small mass external acceleration
end

r = X(1:3,1); % SC position
v = X(4:6,1); % SC velocity

v_dot = [0,2,0;-2,0,0;0,0,0]*v+CR3BP_dUdr(r,mu)+U;

rv_dot = [v;v_dot];

X_dot = rv_dot;
end