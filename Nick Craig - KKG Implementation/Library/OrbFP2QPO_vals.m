function [fv, A] = OrbFP2QPO_vals(X0, rho, eigvecs, eigvals, Nc, u1, u2, w)

X0 = X0(:);
u1 = u1(:);
u2 = u2(:);
eigvals = eigvals(:);


% Order Eigenvalues/vectors in ascending angle and weight
lams = angle(eigvals);

Q_zeros = zeros(6*(2*Nc+1),1);
a0 = sqrt(2) * X0;
Q0 = Q_zeros;
Q0(1:6) = a0; % Create Q0 matrix

% Choose / blend eigenpair
if isscalar(eigvals)
    % Grab and store eigenvector and rotation term
    a1 = real(eigvecs(:,1));
    b1 = imag(eigvecs(:,1));
    tau = lams(1);
    
    % Create vector cooresponding to first frquency
    ab_idx = 1;
    Q1 = Q_zeros;
    idx_a = 6 + 12*(ab_idx-1) + (1:6);
    idx_b = 6 + 12*(ab_idx-1) + (7:12);
    Q1(idx_a,1) = rho*a1; % Set a1
    Q1(idx_b,1) = rho*b1; % Set b1
    
    Q = Q0+Q1;
else
    frequency_ratio = lams(2)/lams(1);
    nm_vals = 1:ceil(Nc/2);
    rationals_mtrx = (1./nm_vals')*nm_vals;
    close_mtrx = abs(rationals_mtrx-frequency_ratio);
    [m,n] = find(close_mtrx == min(min(close_mtrx)),1)

    V1 = eigvecs(:,1);
    V2 = eigvecs(:,2);

    am = (w-1)*real(V1);
    bm = (w-1)*imag(V1);

    an = w*real(V2);
    bn = w*imag(V2);

    % Set a_m and b_m
    Qm = Q_zeros;
    ab_idx = m;
    idx_a = 6 + 12*(ab_idx-1) + (1:6);
    idx_b = 6 + 12*(ab_idx-1) + (7:12);
    Qm(idx_a,1) = rho*am;
    Qm(idx_b,1) = rho*bm;

    % Set a_n and b_n
    Qn = Q_zeros;
    ab_idx = n;
    idx_a = 6 + 12*(ab_idx-1) + (1:6);
    idx_b = 6 + 12*(ab_idx-1) + (7:12);
    Qn(idx_a,1) = rho*an;
    Qn(idx_b,1) = rho*bn;
    
    Q = Q0+Qm+Qn;
    tau = lams(1)/m;
end

% Ok we will pick this up in the morning, but you are hungry. Currenty we
% are trying to figure out how to include terms in the fourier expansion
% that include a secondary principle frequency. This would occur in the
% case where the center space if 2 dimensional. Currently thinking through
% a way to inform the choice in higher order terms in the current
% formulation.

% Free variables
fv = [Q; tau];

% Projected signed area
M = Fourier_AreaMtrx(u1, u2, Nc);
A = Q' * M * Q;

end