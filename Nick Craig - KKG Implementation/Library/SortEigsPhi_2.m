function [vec, val] = SortEigsPhi_2(A, prev_vecs, opts)
% This funciton performs an eigenvalue decomposition on the monodromy
% matrix of a simply periodic orbit. If nothing else if provided, this
% function sorts the eigenvalues in the following order (Reciprocal pairs
% given next to each other).
% 1. Vals and Vecs paired as reciprocals
% 2. Order by |λ|≈1 first, then angle off real.
% 3. Order by real vals next, then complex.
% If a matrix of previous vectors is provided, then this function will
% order the eigenvalues/vectors such that maximum alignment with the
% previous ones are maintained.

if ~exist('prev_vecs','var'); prev_vecs = []; end
if ~exist('opts','var'); opts = []; end
if ~isfield(opts,'tol_one'); opts.tol_one = 1e-4; end
if ~isfield(opts,'tol_unit'); opts.tol_unit = 1e-4; end
if ~isfield(opts,'tol_recip'); opts.tol_recip = 1e-8; end

tol_one = opts.tol_one;
tol_unit = opts.tol_unit;
tol_recip = opts.tol_recip;

[V,D] = eig(A);
w = diag(D);
n = length(w);

% Fix phases
for j = 1:n
    [~,idx] = max(abs(V(:,j)));
    if abs(V(idx,j)) > 0
        V(:,j) = V(:,j)*exp(-1i*angle(V(idx,j)));
    end
end

if isempty(prev_vecs)

    used = false(n,1);
    order = zeros(n,1);
    kk = 1;

    % ---- Special-case pair: zero with one ----
    idx_zero = find(abs(w) < tol_recip & ~used);
    idx_one  = find(abs(w - 1) < tol_one & ~used);

    if ~isempty(idx_zero) && ~isempty(idx_one)
        i0 = idx_zero(1);
        i1 = idx_one(1);

        order(kk:kk+1) = [i1; i0];   % put 1 first, 0 second
        used(i1) = true;
        used(i0) = true;
        kk = kk + 2;
    end

    % ---- General reciprocal pairing ----
    for i = 1:n
        if used(i), continue; end

        lam = w(i);
        cand = find(~used);
        cand(cand == i) = [];

        if isempty(cand)
            order(kk) = i;
            kk = kk + 1;
            used(i) = true;
            break
        end

        [~,loc] = min(abs(w(cand) - 1/lam));
        j = cand(loc);

        order(kk:kk+1) = [i; j];
        used(i) = true;
        used(j) = true;
        kk = kk + 2;
    end

    % Trim in case odd leftover
    order = order(order > 0);

    w = w(order);
    V = V(:,order);

    % ---- Order pairs ----
    lead = w(1:2:end);
    pair_inds = (1:length(lead)).';

    is_one  = (abs(lead-1) < tol_one) | (abs(lead+1) < tol_one);
    is_unit = (abs(abs(lead)-1) < tol_unit) & ~is_one;
    is_real = (abs(imag(lead)) < tol_unit) & ~is_one;
    is_other = ~(is_one | is_unit | is_real);

    pair_order = [pair_inds(is_one); ...
                  pair_inds(is_unit); ...
                  pair_inds(is_real); ...
                  pair_inds(is_other)];

    new_order = zeros(length(w),1);
    for k = 1:length(pair_order)
        p = pair_order(k);
        if 2*p <= length(w)
            new_order(2*k-1:2*k) = [2*p-1; 2*p];
        else
            new_order(2*k-1) = 2*p-1;
        end
    end
    new_order = new_order(new_order > 0);

    val = w(new_order);
    vec = V(:,new_order);

else
    P = prev_vecs ./ max(vecnorm(prev_vecs,2,1), 1e-30);
    Q = V         ./ max(vecnorm(V,2,1),         1e-30);

    S = abs(P' * Q);

    order = zeros(1,n);
    used = false(1,n);

    for i = 1:n
        s = S(i,:);
        s(used) = -inf;
        [~,j] = max(s);
        order(i) = j;
        used(j) = true;
    end

    val = w(order);
    vec = V(:,order);

    for j = 1:n
        alpha = prev_vecs(:,j)' * vec(:,j);
        if abs(alpha) > 0
            vec(:,j) = vec(:,j)*exp(-1i*angle(alpha));
        end
    end
end

end