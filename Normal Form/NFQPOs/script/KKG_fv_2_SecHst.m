function X_sec_hst = KKG_fv_2_SecHst(fv)

fv = nonzeros(fv);
Q = fv(1:end-1,1);
Nh = (length(Q)/6-1)/2;

N_pts = 1000;

th_hst = linspace(0,2*pi,N_pts);
X_sec_hst = zeros(6,N_pts);

for i = 1:N_pts
    X_sec_hst(:,i) = Fourier_Q2X(th_hst(i),Nh)*Q;
end
end

function A = Fourier_Q2X(theta,Nc)

I6 = eye(6);
A = 1/sqrt(2)*I6;

for n = 1:Nc
    A = [A,cos(n*theta)*I6,sin(n*theta)*I6];
end

end
