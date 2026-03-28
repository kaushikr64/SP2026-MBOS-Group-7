function A = Fourier_Q2X(theta,Nc)

I6 = eye(6);
A = 1/sqrt(2)*I6;

for n = 1:Nc
    A = [A,cos(n*theta)*I6,sin(n*theta)*I6];
end

end