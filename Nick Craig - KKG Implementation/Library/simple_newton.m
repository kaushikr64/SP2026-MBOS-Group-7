function [X, info] = simple_newton(fun, X0, opts)
%SIMPLE_NEWTON Minimal Newton solver
%
% Inputs:
%   fun  : function handle -> [F, J] = fun(X)
%   X0   : initial guess
%   opts : struct with fields:
%       max_iter   (default 20)
%       tol        (default 1e-8)
%       step_tol   (default 1e-10)
%       max_step   (default inf)
%
% Outputs:
%   X     : solution
%   info  : struct with convergence info

% ---- defaults ----
if ~exist('opts','var'), opts = struct; end
if ~isfield(opts,'max_iter'), opts.max_iter = 20; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'step_tol'), opts.step_tol = 1e-14; end
if ~isfield(opts,'max_step'), opts.max_step = inf; end

X = X0;

info.converged = false;
info.iter = 0;

for k = 1:opts.max_iter
    
    [F, J] = fun(X);
    
    [res,idx] = max(abs(F)-opts.tol);
    
    fprintf('Iter %d: max|F| = %.3e, idx max = %.0i\n', k-1, res,idx);
    
    % ---- convergence check ----
    if res < 0
        info.converged = true;
        break
    end
    
    % ---- Newton step ----
    dX = -lsqminnorm(J,F);
    % dX = -J\F;
    % dX = -J' * ((J*J') \ F);
    
    % ---- step limiting ----
    step_norm = norm(dX);
    if step_norm > opts.max_step
        dX = dX * (opts.max_step / step_norm);
    end
    
    % ---- step size check ----
    if norm(dX) < opts.step_tol
        fprintf('Step tolerance reached.\n');
        break
    end
    
    % ---- update ----
    X = X + dX;
    
    info.iter = k;
end

end