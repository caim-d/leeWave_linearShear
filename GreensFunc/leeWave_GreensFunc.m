set(0, 'DefaultLineLineWidth', 1.5);
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%
correction = 'linear';

H = 3000;
h0 = 5;
sigma = 250;
f = -1.2e-4;
N = 1e-3;
U0 = 0.1;
UH = 0.1;
Lam = (UH - U0)/H;
LamT = UH - U0;

xStart = -2000;
xEnd = 8000;

num_x = 2^9;
x = linspace(xStart,xEnd,num_x);
dx = x(2)-x(1);
num_z = 2^10;

nMax = 30;

% vector of vertical wavenumbers
m_ = (1:nMax)'*pi/H;
mnMax = m_(end);
dz = H/(num_z - 1);
while 2*pi/mnMax < 10*dz    % make sure vertical wavelengths are well resolved for plotting
    num_z = num_z*2;
    dz = H/(num_z - 1);
end
z = linspace(0,H,num_z);

syms xs U0s sigmas h0s
% Gaussian derivative
hs = -h0s*xs/sigmas*exp(-xs^2/(2*sigmas^2)) * exp(1/2);    % scaled so maximum absolute value is h0
kstar = pi/(2*sigma);

w0s = U0s*diff(hs,xs);
w0_xs = diff(w0s,xs);
w0_x2s = diff(w0_xs,xs);
w0_x3s = diff(w0_x2s,xs);
w0_x4s = diff(w0_x3s,xs);

w0_Ints = U0s*hs;
w0_Int2s = int(w0_Ints,xs);
w0_Int3s = int(w0_Int2s,xs);

h = double(subs(hs,{h0s,sigmas,xs},{h0,sigma,x}));
w0 = double(subs(w0s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_x = double(subs(w0_xs,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_x2 = double(subs(w0_x2s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_x3 = double(subs(w0_x3s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_x4 = double(subs(w0_x4s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));

w0_Int = double(subs(w0_Ints,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_Int2 = double(subs(w0_Int2s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));

w0_Int3 = double(subs(w0_Int3s,{h0s,U0s,sigmas,xs},{h0,U0,sigma,x}));
w0_Int3 = w0_Int3 - w0_Int3(1);

clear xs U0s sigmas h0s hs w0s w0_xs w0_x2s w0_x3s w0_x4s w0_Ints w0_Int2s w0_Int3s;

% cutoff tolerance for determining support of topography
tol = 1e-10;

for i = 1:num_x
    if abs(h(i)) > tol
        x0max = x(i);
        imax = i;
    end
    if abs(h(num_x - i + 1)) > tol
        x0min = x(num_x - i + 1);
        imin = num_x - i + 1;
    end
end

x0 = reshape(x(imin:imax),[1,1,imax - imin + 1]);
num_x0 = imax - imin + 1;

x0inds = imin:imax;

%% get matrix M and vector R, and organize eigenvalues
[M,M4invRw] = get_M_R(nMax,Lam,LamT,H,U0,N,f,w0,w0_x2,w0_x4);

sizeM = length(M);

% eigenvalues and eigenvectors, M*J = J*diag(d_)
% d_ is vector of eigenvalues
[J,d_] = eig(M,'vector');

% make sure complex eigenvalues that should be negative aren't positive
for i = 1:sizeM
    if real(d_(i)) < 1e-13 && abs(imag(d_(i))) > 1e-10 && real(d_(i)) > 0
        d_(i) = 1i*imag(d_(i));
    end
end

% sort eigenvalues end eigenvectors by real part
[d_,inds] = sort(d_,'ComparisonMethod','real');
J = J(:,inds);

Jinv = inv(J);

% check if any of the eigenvalues are repeated
[mindistance, ind_min] = min(abs(diff(d_)));
disp("Minimum distance between eigenvalues = " + mindistance ...
    + ", log10(cond(JMT)) = " + log10(cond(J)));
if mindistance < 1e-11
    disp("Two eigenvalues are (almost) equal");
    pause;
end

dn_ = d_(d_ <= 0);
dp_ = d_(d_ > 0);

%%
disp("x_max = " + x(end));
tic;

[W_,W_Int3_] = get_wModes(dn_,dp_,x0,J,Jinv,M4invRw(:,x0inds),sizeM,nMax,x,num_x,num_x0);

disp("Norm of difference in w_ between W_ and W_Int3_ = " + norm(W_(1:nMax) - W_Int3_(3*nMax+1:4*nMax)));
disp("Norm of w_ in W_ = " + norm(W_(1:nMax)));

runtime = toc;

% save("UH" + UH + "_tol" + tol + "_num_x" + num_x + "_nMax" + nMax + "_L" + L + ".mat");

%% get wave fields from modes
[xplot,zplot] = meshgrid(x,z);

[w,u,v,b,phi,...
    xi,eta,zeta,zeta_z,...
    w_x,u_x,v_x,b_x,...
    w_z,u_z,v_z,b_z,phi_z,...
    w_x2,u_x2,v_x2,b_x2,...
    w_z2,u_z2,v_z2,b_z2,...
    w_xz,u_xz,v_xz,b_xz,...
    w_x3,u_x3,v_x3,b_x3,...
    w_x2z,u_x2z,v_x2z,b_x2z,...
    w_xz2,u_xz2,v_xz2,b_xz2,...
    w_x4,u_x4,v_x4,b_x4,...
    w_x3z,u_x3z,v_x3z,b_x3z,...
    w_x2z2,u_x2z2,v_x2z2,b_x2z2] = get_waveFields(W_,W_Int3_,m_,w0,w0_x,w0_x2,w0_x3,...
                    w0_Int,w0_Int2,w0_Int3,zplot,nMax,f,N,U0,H,Lam);

xi_x = (Lam*zeta + u)./(U0 + Lam*zplot);
eta_x = v./(U0 + Lam*zplot);
zeta_x = w./(U0 + Lam*zplot);

% linear PV
q = v_x + f/N^2*b_z - f*Lam/N^2*(u_z - w_x);

% get wave properties
% Lagrangian pseudomomentum
pmomL1 = -xi_x.*(u + Lam*zeta - f/2*eta) - eta_x.*(v + f/2*xi) - zeta_x.*w;

% Lagrangian pseudomomentum fluxes
FL1 = (U0 + Lam*zplot).*pmomL1 + 1/2*(u.^2 + v.^2 + w.^2 - N^2*zeta.^2) ...
                               + f/2*(Lam*eta.*zeta + v.*xi - u.*eta) ...
                               + Lam*u.*zeta - xi_x.*phi + Lam^2/2*zeta.^2;
FL3 = - zeta_x.*phi;

% Eulerian pseudomomentum
pmomE1 = 1/N^2*b.*(u_z - w_x);

% Eulerian pseudomomentum fluxes
FE1 = (U0 + Lam*zplot).*pmomE1 + 1/2*(b.^2/N^2 + u.^2 - v.^2 - w.^2);
FE3 = u.*w - f*v.*b/N^2;

% wave energy density
E = 1/2*(u.^2 + v.^2 + w.^2 + b.^2/N^2);

% linear PV
q = v_x + f/N^2*b_z - f*Lam/N^2*(u_z - w_x);

Fr = sqrt(u_z.^2 + v_z.^2)/N;

%%
figure;

xCutoff = 8000;
plotInds = 1:length(x(x <= xCutoff));

field_to_plot = w;

PLOT = pcolor(xplot(:,plotInds), zplot(:,plotInds), field_to_plot(:,plotInds));
set(PLOT, 'EdgeColor', 'none');
title("\(w^\prime\) for \(\Lambda = \) " + Lam + ", \(n_{\mathrm{max}}\) = " + nMax, "\(\sigma = \) " + sigma);
colorbar;
xlabel('\(x\)');
ylabel('\(z\)');
maxclim = max(abs(clim));
clim([-maxclim,maxclim]);
hold on
yline((abs(f)/kstar - U0)/Lam,'--k');   % critical levels
yline((-abs(f)/kstar - U0)/Lam,'--k');
yline((N/kstar - U0)/Lam,'-.k');        % turning level
yline(-U0/Lam,'--g');
hold off

xrange = xlim;
pbaspect([xrange(2)-xrange(1),H,1]);
