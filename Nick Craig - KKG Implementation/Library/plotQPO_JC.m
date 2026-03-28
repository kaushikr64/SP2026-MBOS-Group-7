function fig = plotQPO_JC(Q, Nc, Sys)

if isvector(Q)
    Q = Q(:);
end

hold on;
grid on;

N_orbs = size(Q,2);

for k = 1:N_orbs
    Qk = Q(:,k);
    fplot(@(th) QPO_JCFcn(Qk, Nc, Sys.mu, th), [0, 2*pi]);
end

end

function val = QPO_JCFcn(Q, Nc, mu, theta)

val = arrayfun(@(th) eval_JC(th, Q, Nc, mu), theta);

end

function JC = eval_JC(theta, Q, Nc, mu)

X = Fourier_Q2X(theta, Nc) * Q;
JC = CR3BP_JC(X, mu);

end