function fig = plotQPOFamilyPoinSection(Q, Np, Nc, Sys,section)
ax = gca;
hold on;
grid on;
% axis equal;
N_orbs = size(Q,2);

for k = 1:N_orbs
    Qk = Q(:,k);

    % Continuous Fourier curve in position space
    fplot3(@(th) QPO_PosFcn(Qk,Nc,1,th)*Sys.Ls, ...
           @(th) QPO_PosFcn(Qk,Nc,2,th)*Sys.Ls, ...
           @(th) QPO_PosFcn(Qk,Nc,3,th)*Sys.Ls, ...
           [0, 2*pi]);

    % Optional discrete points used in shooting
    if Np > 0
        Xp = zeros(3,Np);
        for i = 1:Np
            th = (i-1)*2*pi/Np;
            Xth = Fourier_Q2X(th,Nc) * Qk;
            Xp(:,i) = Xth(1:3) * Sys.Ls;
        end
        plot3(Xp(1,:), Xp(2,:), Xp(3,:), 'o', 'MarkerSize', 5);
    end
end

% Orient camera to look along the position part of the section normal
n_pos = section.n_sec(1:3);
if norm(n_pos) > 1e-12
    n_pos = n_pos / norm(n_pos);
    az = atan2d(n_pos(2), n_pos(1));
    el = atan2d(n_pos(3), sqrt(n_pos(1)^2 + n_pos(2)^2));
    view(ax,az+90, el);
end

end

function val = QPO_PosFcn(Q, Nc, idx, theta)

a0 = Q(1:6);
A = reshape(Q(7:end), 6, []);

X = (1/sqrt(2)) * a0;

col = 1;
for n = 1:Nc
    a_n = A(:,col);   col = col + 1;
    b_n = A(:,col);   col = col + 1;
    X = X + a_n*cos(n*theta) + b_n*sin(n*theta);
end

val = X(idx);

end