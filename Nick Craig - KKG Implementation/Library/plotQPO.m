function plotQPO(Q,propagate,Np,Nc,Sys)
for i = 1:Np
    th = i*2*pi/Np;
    A_th0 = Fourier_Q2X(th,Nc);
    X_th0 = A_th0*Q;

    [~,~,~,sol] = propagate(X_th0);

    t_hst = linspace(0,sol.xe(1),2000);
    X_hst = deval(sol,t_hst);

    plot3(X_hst(1,:)*Sys.Ls,X_hst(2,:)*Sys.Ls,X_hst(3,:)*Sys.Ls,'color','b');
    % plot3(X_hst(1,1),X_hst(2,1),X_hst(3,1),'color','b',...
    %     'Marker','o','LineStyle','none');
    % plot3(X_hst(1,end),X_hst(2,end),X_hst(3,end),'color','r',...
    %     'Marker','x','LineStyle','none');
    % fplot3(QPO_PosFcn(Q,Nc,1), QPO_PosFcn(Q,Nc,2), QPO_PosFcn(Q,Nc,3));
    % view(0,0);
    
    hold on;
    grid on;
end
end

function f = QPO_PosFcn(Q, Nc, idx)

a0 = Q(1:6);

A = reshape(Q(7:end), 6, []); % [a1 b1 a2 b2 ...]

f = @(theta) arrayfun(@(th) ...
    eval_theta(th, a0, A, Nc, idx), theta);

end

function val = eval_theta(theta, a0, A, Nc, idx)

X = (1/sqrt(2)) * a0;

col = 1;
for n = 1:Nc
    a_n = A(:,col);   col = col + 1;
    b_n = A(:,col);   col = col + 1;

    X = X + a_n*cos(n*theta) + b_n*sin(n*theta);
end

val = X(idx);

end