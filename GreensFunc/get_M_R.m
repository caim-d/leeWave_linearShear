function [M,M4invRw] = get_M_R(nMax,Lam,LamT,H,U0,N,f,w0,w0_x2,w0_x4)

[A0,A1s,A2s,A3s,~] = get_AMatrices(nMax);
[Bsc,~] = get_BMatrices(nMax);
[s0,s1,s2,s3,s4,~,~] = get_sVectors(nMax);

% vector of m_n
m_ = (1:nMax)'*pi/H;
mstar_ = m_;
for n = 1:2:nMax
    mstar_(n) = -mstar_(n);
end

% assume not hydrostatic and no damping. should work for shear and no shear
M4 = U0^3*A0 + 3*U0^2*LamT*A1s + 3*U0*LamT^2*A2s + LamT^3*A3s;
M2 = N^2*(U0*A0 + LamT*A1s) - M4*diag(m_.^2);
M0 = - f^2*(U0*A0 + LamT*A1s)*diag(m_.^2) - 2*f^2*LamT/H*Bsc*diag(m_);


M = zeros(4*nMax);    % 4*nMax vertical modes
M(1:3*nMax,nMax+1:4*nMax) = eye(3*nMax);
M(3*nMax+1:4*nMax,1:nMax) = -M4\M0;
M(3*nMax+1:4*nMax,2*nMax+1:3*nMax) = -M4\M2;



% find right hand side vectors
R1_4 = U0^3*(s1 - s0) + 3*U0^2*LamT*(s2 - s1) ...
    + 3*U0*LamT^2*(s3 - s2) + LamT^3*(s4 - s3);
R1_2 = N^2*( U0*(s1 - s0) + LamT*(s2 - s1) );
R1_0 = -2*f^2*Lam/H*s0;

R1 = R1_4.*w0_x4 + R1_2.*w0_x2 + R1_0.*w0;

M4invRw = M4\R1;


disp("cond(M_4) = " + cond(M4));

