function W_Int3_j = get_W_Int3_j(dn_,dp_,X,x0,J,Jinv,M4invRw,sizeM,nMax)

num_dn = length(dn_);
num_dp = length(dp_);

% x0 is in dimension 3
x0leX = x0(x0 <= X);
x0geX = x0(x0 >= X);

num_x0leX = length(x0leX);
num_x0geX = length(x0geX);

% the d_ vector is ordered by real part
dn_maxInd = num_dn;
dp_minInd = sizeM - num_dp + 1;

Jn = J(:,1:dn_maxInd);
Jp = J(:,dp_minInd:end);

% -------------------------------------------------------------------------
% first do x0 < X
% contribution from dn_ to x0 < X
Ex0leX_Int3n = 1./((dn_.^3).') .* (exp( dn_.' .* (X - x0leX) ) - 1) ...
            - 1./((dn_.^2).') .* (X - x0leX) ...
            - 1./(2*dn_.') .* ((X - x0leX).^2);
Jn_Ex0leX_Int3n = Jn .* Ex0leX_Int3n;

% contribution from dp_ to x0 < X
Ex0leX_Int3p = -1./((dp_.^3).').*ones(1,1,num_x0leX) - 1./((dp_.^2).').*(X - x0leX) ...
            - 1./(2*dp_.') .* ((X - x0leX).^2);
Jp_Ex0leX_Int3p = Jp .* Ex0leX_Int3p;

Gx0leX_M4invRw = pagemtimes([Jn_Ex0leX_Int3n,Jp_Ex0leX_Int3p],Jinv(:,3*nMax+1:4*nMax));
Gx0leX_M4invRw = real(Gx0leX_M4invRw);

% -------------------------------------------------------------------------
% now do x0 > X
Ex0geX_Int3 = -1./((dp_.^3).') .* exp( dp_.' .* (X - x0geX) );
Jp_Ex0geX_Int3 = Jp .* Ex0geX_Int3;

Gx0geX_M4invRw = pagemtimes(Jp_Ex0geX_Int3,Jinv(dp_minInd:end,3*nMax+1:4*nMax));
Gx0geX_M4invRw = real(Gx0geX_M4invRw);

% -------------------------------------------------------------------------
integrand_x0leX = pagemtimes(Gx0leX_M4invRw,M4invRw(:,:,x0 <= X));
integrand_x0geX = pagemtimes(Gx0geX_M4invRw,M4invRw(:,:,x0 >= X));

W_Int3_j = trapz(reshape(x0leX,[num_x0leX,1]),integrand_x0leX,3) ...
        + trapz(reshape(x0geX,[num_x0geX,1]),integrand_x0geX,3);
