function fig = plotQPOOrbit_3D(Q, propagate, Np, Nc, Sys)

fig = figure();
ax = axes(fig);
hold on;
grid on;
axis equal;
view(ax,-20,20);

for i = 1:Np
    th = (i-1)*2*pi/Np;
    X0 = Fourier_Q2X(th, Nc) * Q;

    [~,~,~,sol] = propagate(X0);

    t_hst = linspace(0, sol.xe(1), 2000);
    X_hst = deval(sol, t_hst);

    plot3(Sys.Ls * X_hst(1,:), ...
          Sys.Ls * X_hst(2,:), ...
          Sys.Ls * X_hst(3,:), ...
          'b');
end

end