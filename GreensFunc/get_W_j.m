function W_j = get_W_j(dn_,dp_,X,x0,J,Jinv,M4invRw,sizeM,nMax)

num_dn = length(dn_);   % number of eigenvalues which are nonpositive
num_dp = length(dp_);   % and positive

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
Ex0leX = exp( dn_.' .* (X - x0leX) ); % dn_.' here so that this vector is compatible dimension to multiply Jn
Jn_Ex0leX = Jn .* Ex0leX;

Gx0leX_M4invRw = pagemtimes(Jn_Ex0leX,Jinv(1:dn_maxInd,3*nMax+1:4*nMax));
Gx0leX_M4invRw = real(Gx0leX_M4invRw);

% -------------------------------------------------------------------------
% now do x0 > X
Ex0geX = -exp( dp_.' .* (X - x0geX) );   % this one has minus in front
Jp_Ex0geX = Jp .* Ex0geX;

Gx0geX_M4invRw = pagemtimes(Jp_Ex0geX,Jinv(dp_minInd:end,3*nMax+1:4*nMax));
Gx0geX_M4invRw = real(Gx0geX_M4invRw);

% -------------------------------------------------------------------------
integrand_x0leX = pagemtimes(Gx0leX_M4invRw,M4invRw(:,:,x0 <= X));
integrand_x0geX = pagemtimes(Gx0geX_M4invRw,M4invRw(:,:,x0 >= X));

W_j = trapz(reshape(x0leX,[num_x0leX,1]),integrand_x0leX,3) ...
        + trapz(reshape(x0geX,[num_x0geX,1]),integrand_x0geX,3);
