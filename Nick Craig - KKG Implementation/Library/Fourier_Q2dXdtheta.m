function B = Fourier_Q2dXdtheta(theta,Nc)

I6 = eye(6);
B = zeros(6);

for n = 1:Nc
    B = [B,-n*sin(n*theta)*I6,n*cos(n*theta)*I6];
end

end