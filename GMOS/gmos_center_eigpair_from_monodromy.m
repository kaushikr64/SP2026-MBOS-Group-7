function cen = gmos_center_eigpair_from_monodromy(tspan, xPO0, mu, mode)
%GMOS_CENTER_EIGPAIR_FROM_MONODROMY Select planar or vertical center eigenpair from monodromy.
%
% Inputs
%   tspan : [t0 tf] (one period)
%   xPO0  : 6x1 periodic orbit initial condition
%   mu    : CR3BP mass parameter
%   mode  : 'planar' (default) or 'vertical'
%
% Output
%   cen.lambda : selected eigenvalue (imag(lambda) > 0)
%   cen.w      : eigenvector
%   cen.rho0   : angle(lambda)
%   cen.idx    : selected index in eig output
%   cen.proj   : planar/vertical projection norms for the selected vector

    if nargin < 4 || isempty(mode)
        mode = 'planar';
    end

    M_eigs = monodromy_CR3BP(tspan, xPO0, mu);
    V = M_eigs(:, 1:6);
    D = M_eigs(:, 7:12);
    lam = diag(D);

    tol_unit = 1e-6;
    tol_imag = 1e-10;

    % Candidates near the unit circle with nontrivial imaginary part
    cand = find(abs(abs(lam) - 1) < tol_unit & abs(imag(lam)) > tol_imag);

    if isempty(cand)
        error('No complex unit-magnitude eigenvalue candidates found.');
    end

    % Deterministic member of each conjugate pair: imag(lambda) > 0
    cand = cand(imag(lam(cand)) > 0);

    if isempty(cand)
        error('No candidates with imag(lambda) > 0 found.');
    end

    % Compute projection scores for each candidate
    planar_score = zeros(numel(cand), 1);
    vertical_score = zeros(numel(cand), 1);

    for k = 1:numel(cand)
        w = V(:, cand(k));
        planar_score(k) = norm([w(1); w(2); w(4); w(5)]);
        vertical_score(k) = norm([w(3); w(6)]);
    end

    % Select by requested mode
    if strcmpi(mode, 'vertical')
        [bestScore, kbest] = max(vertical_score);
    else
        [bestScore, kbest] = max(planar_score);
    end

    i0 = cand(kbest);

    lambda = lam(i0);
    w = V(:, i0);

    % Package output
    cen = struct();
    cen.lambda = lambda;
    cen.w = w;
    cen.rho0 = angle(lambda);
    cen.idx = i0;

    cen.proj = struct();
    cen.proj.planar = norm([w(1); w(2); w(4); w(5)]);
    cen.proj.vertical = norm([w(3); w(6)]);

    % Warn if requested mode is essentially absent
    if bestScore < 0.2
        warning('Requested mode "%s" is weak (score %.3e). You may not have a meaningful %s center mode here.', ...
            mode, bestScore, mode);
    end
end