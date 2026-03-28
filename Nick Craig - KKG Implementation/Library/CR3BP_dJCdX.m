function dJC = CR3BP_dJCdX(X,mu)
r = X(1:3,1);
v = X(4:6,1);

dJC = [2*CR3BP_dUdr(r,mu)',-2*v'];
end