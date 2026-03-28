function h = PoincareSec(X, X_sec, n_sec)
% Definition of a general Poincare section and when a state lies on it
% (zero).
X_sec = X_sec(:);
n_sec = n_sec(:);

h = dot(n_sec*ones(size(X(1,:))),X-X_sec);
end