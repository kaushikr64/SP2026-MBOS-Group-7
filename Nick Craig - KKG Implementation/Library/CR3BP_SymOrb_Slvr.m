function [X_cnv,DFDX_cnv] = CR3BP_SymOrb_Slvr(X_g,mu,const_var)
% The purpose of this function is to solve for an initial state in the
% CR3BP system that results in a XZ-plane symmetric periodic orbit. The
% function uses a "good" initial guess for the x,z,vy,T/2 and attemepts to
% solve the non-linear shooting function nearest to the initial guess.
% Given the nature of the symmetric orbits, one could choose to hold either
% the initial x or z position guess and solve for the other values. This is
% done through the second inpupt. The outputs of the function consist of
% the converged non-zero initial state values and half period, along with
% the jacobian of the shooting function at the final time with respect to
% the non-zero state and half period values. This jacobian is useful when
% propagating a family of symmetric orbits.
% 
% Inputs:
%   X_g [x,z,ydot,T_2]'
%       x:      initial guess x position
%       z:      ...           z position
%       ydot:   ...           y velocity
%       T2:     ...           half period
%   mu:         gravitational constant of the 2 system primaries
%   const_var:  "x","z" or NaN specifying the coordinate to hold constant
% 
% Outputs:
%   X_cnv [x_cnv,z_cnv,ydot_cnv,T_2_cnv]'
%       x_cnv:      initial converged x position
%       z_cnv:      ...               z position
%       ydot_cnv:   ...               y velocity
%       T_2_cnv:    ...               half period
%  DFDX_cnv
% 
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

if ~exist('mu','var')
    mu = 1.215059e-2; % Earth-Moon MU
end

int_tol = 2.3e-14; % Integration tolerance used inside of shooting function
sf_tol = 3e-14; % Shooting funcion tolerance
fslveopts = optimoptions('fsolve','FunctionTolerance',sf_tol,'StepTolerance',sf_tol,...
    'OptimalityTolerance',sf_tol,'SpecifyObjectiveGradient',true,'Algorithm',...
    'levenberg-marquardt','Display','none'); % fsolve options

if ~exist('const_var','var')
    Z_g = X_g; % Intial free-variable guess

    Z_cnv = fsolve(@(Z) PrdSymmOrb_SF_Gen(Z,mu),Z_g,fslveopts);
    [~,~,DFDX_cnv] = PrdSymmOrb_SF_Gen(Z_cnv,mu);    
    X_cnv = Z_cnv;
elseif isequal(const_var,"x")
    x = X_g(1);
    Z_g = X_g(2:4);
    
    Z_cnv = fsolve(@(Z) PrdSymmOrb_SF_Xcnst(x,Z,mu),Z_g,fslveopts);
    [~,~,DFDX_cnv] = PrdSymmOrb_SF_Xcnst(x,Z_cnv,mu);
    X_cnv = [x;Z_cnv];    
elseif isequal(const_var,"z")
    z = X_g(2);
    Z_g = [X_g(1);X_g(3:4)];

    Z_cnv = fsolve(@(Z) PrdSymmOrb_SF_Zcnst(z,Z,mu),Z_g,fslveopts);
    [~,~,DFDX_cnv] = PrdSymmOrb_SF_Zcnst(z,Z_cnv,mu);
    X_cnv = [Z_cnv(1);z;Z_cnv(2:3)];
end

function [F,DFDZ,DFDX] = PrdSymmOrb_SF_Gen(Z,mu)
% This is the CR3BP symmetric orbit shooting function. The input is an
% initial guess for the x/z coordinates, the y velocity, and the half
% period given in a column vector. The shooting function itself, propagates
% the initial state through the half period and returns the final y
% position and x/z velocity. For a true symmetric orbit it the CR3BP
% system, these 3 values will be zero. The objective is to iteratively
% solve the shooting function to find an initial guess that results in
% these values being zero.
% 
% This is the general shooting function that solves for all 4 values.

% Extraction of values from free variable
rv_0 = [Z(1),0,Z(2),0,Z(3),0]'; % Initial state from free variables
T2 = Z(4); % Half period to integrate through
opts = odeset('RelTol',int_tol,'AbsTol',int_tol); % Integration options

% Integration of state through half period and extraction of final values
ode = ode45(@(t,X) CR3BP_STM_EOM(X,mu),[0,T2],[rv_0;...
    reshape(eye(6),[36,1])],opts); % Integrated trajectory
X_phi_f = deval(ode,T2); % Final state and STM vector
X_f = X_phi_f(1:6,1); % Final state vector
Xdot_f = CR3BP_EOM(X_f,mu); % Final state derivative
phi_f = reshape(X_phi_f(7:42,1),[6,6]); % Final STM

% Shooting function evaluation (want this to be zero)
F = [X_f(2,1);X_f(4,1);X_f(6,1)];

% Derivative of shooting function with respect to free variables
% [x,z,ydot,T/2]
DFDZ = [phi_f(2,1),phi_f(2,3),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,1),phi_f(4,3),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,1),phi_f(6,3),phi_f(6,5),Xdot_f(6,1)];

% Derivative of the shooting function with respect to the outside free
% variables [x,z,ydot,T/2]
DFDX = [phi_f(2,1),phi_f(2,3),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,1),phi_f(4,3),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,1),phi_f(6,3),phi_f(6,5),Xdot_f(6,1)];
end

function [F,DFDZ,DFDX] = PrdSymmOrb_SF_Xcnst(x,Z,mu)
% This is the CR3BP symmetric orbit shooting function. The input is an
% initial guess for the x/z coordinates, the y velocity, and the half
% period given in a column vector. The shooting function itself, propagates
% the initial state through the half period and returns the final y
% position and x/z velocity. For a true symmetric orbit it the CR3BP
% system, these 3 values will be zero. The objective is to iteratively
% solve the shooting function to find an initial guess that results in
% these values being zero.
% 
% This is the general shooting function that solves for all 4 values.

% Extraction of values from free variable
rv_0 = [x,0,Z(1),0,Z(2),0]'; % Initial state from free variables
T2 = Z(3); % Half period to integrate through
opts = odeset('RelTol',int_tol,'AbsTol',int_tol); % Integration options

% Integration of state through half period and extraction of final values
ode = ode45(@(t,X) CR3BP_STM_EOM(X,mu),[0,T2],[rv_0;...
    reshape(eye(6),[36,1])],opts); % Integrated trajectory
X_phi_f = deval(ode,T2); % Final state and STM vector
X_f = X_phi_f(1:6,1); % Final state vector
Xdot_f = CR3BP_EOM(X_f,mu); % Final state derivative
phi_f = reshape(X_phi_f(7:42,1),[6,6]); % Final STM

% Shooting function evaluation (want this to be zero)
F = [X_f(2,1);X_f(4,1);X_f(6,1)];

% Derivative of shooting function with respect to free variables
% [z,ydot,T/2]
DFDZ = [phi_f(2,3),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,3),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,3),phi_f(6,5),Xdot_f(6,1)];

% Derivative of the shooting function with respect to the outside free
% variables [x,z,ydot,T/2]
DFDX = [phi_f(2,1),phi_f(2,3),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,1),phi_f(4,3),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,1),phi_f(6,3),phi_f(6,5),Xdot_f(6,1)];
end

function [F,DFDZ,DFDX] = PrdSymmOrb_SF_Zcnst(z,Z,mu)
% This is the CR3BP symmetric orbit shooting function. The input is an
% initial guess for the x/z coordinates, the y velocity, and the half
% period given in a column vector. The shooting function itself, propagates
% the initial state through the half period and returns the final y
% position and x/z velocity. For a true symmetric orbit it the CR3BP
% system, these 3 values will be zero. The objective is to iteratively
% solve the shooting function to find an initial guess that results in
% these values being zero.
% 
% This is the general shooting function that solves for all 4 values.

% Extraction of values from free variable
rv_0 = [Z(1),0,z,0,Z(2),0]'; % Initial state from free variables
T2 = Z(3); % Half period to integrate through
opts = odeset('RelTol',int_tol,'AbsTol',int_tol); % Integration options

% Integration of state through half period and extraction of final values
ode = ode45(@(t,X) CR3BP_STM_EOM(X,mu),[0,T2],[rv_0;...
    reshape(eye(6),[36,1])],opts); % Integrated trajectory
X_phi_f = deval(ode,T2); % Final state and STM vector
X_f = X_phi_f(1:6,1); % Final state vector
Xdot_f = CR3BP_EOM(X_f,mu); % Final state derivative
phi_f = reshape(X_phi_f(7:42,1),[6,6]); % Final STM

% Shooting function evaluation (want this to be zero)
F = [X_f(2,1);X_f(4,1);X_f(6,1)];

% Derivative of shooting function with respect to free variables
% [x,ydot,T/2]
DFDZ = [phi_f(2,1),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,1),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,1),phi_f(6,5),Xdot_f(6,1)];

% Derivative of the shooting function with respect to the outside free
% variables [x,z,ydot,T/2]
DFDX = [phi_f(2,1),phi_f(2,3),phi_f(2,5),Xdot_f(2,1);...
    phi_f(4,1),phi_f(4,3),phi_f(4,5),Xdot_f(4,1);...
    phi_f(6,1),phi_f(6,3),phi_f(6,5),Xdot_f(6,1)];
end
end

