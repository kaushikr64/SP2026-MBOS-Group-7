function [F, DF] = CR3BP_QPO_SF(fv_guess, opts)

% System parameters
mu = opts.mu;

% Length terms
Nx = 6;
Np = opts.Np; % Number of points on Poincare section to propagate
Nc = opts.Nc; % Number of coefficients to use in Fourier approximation

% Poincare section values
X_sec = opts.section.X_sec;
n_sec = opts.section.n_sec;
null_sec = null(n_sec');

% Area constaint terms
u1 = opts.con_area.u1;
u2 = opts.con_area.u2;
con_A = opts.con_area.A;

% Jacobi constant constraint terms
con_JC = opts.con_JC.JC;

% Define the free variables
fv_Q = fv_guess(1:Nx*(1+2*Nc),1); % Fourier coefficients
fv_tau = fv_guess(Nx*(1+2*Nc)+1,1); % Angular twist of the toroid on return

% Initialize function value and gradient matrix
F = [];
DF = [];

% Loop through the number of points, propagate, grab derivatives and all

for i = 0:(Np-1)
    % Create the ith point on the Fourier curve
    th = i*2*pi/Np;
    A_th0 = Fourier_Q2X(th,Nc);
    A_th0tau = Fourier_Q2X(th+fv_tau, Nc);
    B_th0tau = Fourier_Q2dXdtheta(th+fv_tau, Nc);

    X_th0 = A_th0*fv_Q;
    X_th0tau = A_th0tau*fv_Q;

    [X_return, DP_return, ~, ~] = CR3BP_Prop2Poincare(X_th0, mu, opts.section,[0,20]);

    F_return = null_sec'*(X_return-X_th0tau);
    F_init = n_sec'*(X_th0 - X_sec);

    F = [F;F_return;F_init];

    dFdQ_return = null_sec'*(DP_return*A_th0 - A_th0tau);
    dFdtau_return = -null_sec'*B_th0tau*fv_Q;
    dFdQ_init = n_sec'*A_th0;
    dFdtau_init = 0;

    DF = [DF;...
        dFdQ_return, dFdtau_return;...
        dFdQ_init, dFdtau_init];
end

% Now for the Jacobi constant terms
F_JC = CR3BP_JC(X_th0,mu) - con_JC;
F = [F;F_JC];
dFdQ_JC = CR3BP_dJCdX(X_th0,mu)*A_th0;
dFdtau_JC = 0;
DF = [DF;...
    dFdQ_JC, dFdtau_JC];

% And lastly, the Area terms
M = Fourier_AreaMtrx(u1, u2, Nc);
scl_trm = abs(con_A);
F_area = (fv_Q'*M*fv_Q - con_A)./scl_trm;
F = [F;F_area];

dFdQ_area = (2/scl_trm)*fv_Q'*M;
dFdtau_area = 0;
DF = [DF;...
    dFdQ_area, dFdtau_area];



end