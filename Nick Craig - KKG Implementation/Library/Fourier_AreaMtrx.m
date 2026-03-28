function M = Fourier_AreaMtrx(u1,u2,Nc)
S = u1*u2' - u2*u1';

blocks{1} = zeros(6);
for n = 1:Nc
    blocks{n+1} = (n*pi/2)*[zeros(6),S;...
        -S, zeros(6)];
end
M = blkdiag(blocks{:});
end