function dfv = sol2nextfv(dFdfv_sol, A_init, A_scale)

grad_partial = dFdfv_sol(1:end-1, :);
grad_A = dFdfv_sol(end, :);

[~, S, V] = svd(grad_partial, 'econ');
s = diag(S);

% Pick near-null subspace
tol = max(size(grad_partial)) * eps(max(s));
idx = find(s <= max(10*tol, 1e-8*max(s)));

% If none are tiny enough, fall back to the smallest few
if isempty(idx)
    k = min(3, size(V,2));
    T = V(:, end-k+1:end);
else
    T = V(:, idx);
end

% Direction in that subspace maximizing area increase
c = (grad_A * T).';
t = T * c;

% If projected area gradient is tiny, fail loudly
if norm(t) < 1e-14
    error('Projected area gradient is near zero; continuation direction is ill-defined.');
end

% Normalize
t = t / norm(t);

% Orient so area increases
if grad_A * t < 0
    t = -t;
end

DeltaA_target = A_init * (A_scale - 1);
den = grad_A * t;

if abs(den) < 1e-12
    error('Directional derivative of area is too small; step would blow up.');
end

alpha = DeltaA_target / den;
dfv = alpha * t;

end