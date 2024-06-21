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
UH = 0.01;
Lam = (UH - U0)/H;
LamT = UH - U0;
kstar = pi/(2*sigma);

nu = 1e-2;

xStart = -2000;
xEnd = 8000;

num_x = 2^9;
x = linspace(xStart,xEnd,num_x+1);
x = x(1:end-1);
dx = x(2)-x(1);
num_z = 2^10;

nMax = 240;

% horizontal wavenumbers
k = [0:num_x/2-1, -num_x/2:-1];
k = k*2*pi/(xEnd - xStart);
kdim3 = reshape(k,[1,1,num_x]);


%%
% vector of vertical wavenumbers
m_ = (1:nMax)'*pi/H;
mnMax = m_(end);
dz = H/(num_z - 1);
while 2*pi/mnMax < 10*dz    % make sure vertical wavelengths are well resolved for plotting
    num_z = num_z*2;
    dz = H/(num_z - 1);
end
z = linspace(0,H,num_z);

% mean 0 topography
h = -h0*x/sigma.*exp(-x.^2/(2*sigma^2)) * exp(1/2);
hk = fft(h);

w0k = U0*1i*k.*hk;
w0_Intk = U0*hk;
u0k = 1/H*w0_Intk;

w0 = real(ifft(w0k));
w0_Int = real(ifft(w0_Intk));
u0 = real(ifft(u0k));

w0kdim3 = reshape(w0k,[1,1,num_x]);
u0kdim3 = reshape(u0k,[1,1,num_x]);

%% get matrix M and vector R, and organize eigenvalues
[M,R] = get_M_R(k,nMax,Lam,LamT,H,U0,N,f,nu,w0k);

% get wave modes
Wk_ = pagemldivide(M,R);

W_Intk_ = zeros(nMax,1,num_x);
W_Intk_(:,:,2:end) = -1i./kdim3(2:end).*Wk_(:,:,2:end);

Uk_ = - m_.*W_Intk_;

W_ = real(ifft(Wk_,num_x,3));

U_ = real(ifft(Uk_,num_x,3));

% get V_
Mv = zeros(nMax+1,nMax+1,num_x);
% get needed matrices and vectors
[A0,A1s,~,~,A1c] = get_AMatrices(nMax);
[s0,s1,s2,~,~,~,~] = get_sVectors(nMax);
c1 = get_cVectors(nMax);

Mv(1,1,:) = 1i*kdim3 * (U0 + LamT/2);
Mv(1,2:end,:) = 1i*LamT/2 * kdim3 .* c1';
Mv(2:end,1,:) = 1i*LamT * kdim3 .* c1;
Mv(2:end,2:end,:) = 1i*kdim3 .* (U0*A0 + LamT*A1c) + nu*diag(m_.^2);

Vk_ = zeros(nMax+1,1,num_x);
Vk_(:,:,2:end) = pagemldivide(Mv(:,:,2:end),-f*cat(1,u0kdim3(:,:,2:end),Uk_(:,:,2:end)));
V_ = real(ifft(Vk_,num_x,3));

% get B_
% first Bc_
Mbc = zeros(nMax+1,nMax+1,num_x);
% get needed matrices and vectors

Mbc(1,1,:) = 1i*kdim3 * (U0 + LamT/2);
Mbc(1,2:end,:) = 1i*LamT/2 * kdim3 .* c1';
Mbc(2:end,1,:) = 1i*LamT * kdim3 .* c1;
Mbc(2:end,2:end,:) = 1i*kdim3 .* (U0*A0 + LamT*A1c) + nu*diag(m_.^2);

Bck_ = zeros(nMax+1,1,num_x);
Bck_(:,:,2:end) = pagemldivide(Mbc(:,:,2:end),f*Lam*Vk_(:,:,2:end));
Bc_ = real(ifft(Bck_,num_x,3));

%second Bs_
Mbs = zeros(nMax+1,nMax+1,num_x);

Mbs(1,1,:) = 1i*kdim3 * (U0 + LamT/3);
Mbs(1,2:end,:) = 1i* kdim3 .* (U0*s0' + LamT*s1') + nu*(s0'*diag(m_.^2));
Mbs(2:end,1,:) = 1i* kdim3 .* (U0*(s0 - s1) + LamT*(s1 - s2));
Mbs(2:end,2:end,:) = 1i*kdim3 .* (U0*A0 + LamT*A1s) + nu*diag(m_.^2);

Rbs = -N^2*cat(1,w0kdim3 + pagemtimes(s0',Wk_),Wk_ + (s0 - s1).*w0kdim3);

Bsk_ = zeros(nMax+1,1,num_x);
Bsk_(:,:,2:end) = pagemldivide(Mbs(:,:,2:end),Rbs(:,:,2:end));
Bs_ = real(ifft(Bsk_,num_x,3));

% reshape all here
W_ = reshape(W_,[nMax,num_x]);
U_ = reshape(U_,[nMax,num_x]);
V_ = reshape(V_,[nMax+1,num_x]);
Bc_ = reshape(Bc_,[nMax+1,num_x]);
Bs_ = reshape(Bs_,[nMax+1,num_x]);


%% get wave fields from modes
[xplot,zplot] = meshgrid(x,z);

w = 0*zplot;
u = w;
v = w;
b = w;

for n = 1:nMax
    mn = m_(n);

    wn = W_(n,:);

    un = U_(n,:);
    vn = V_(n+1,:);
    bcn = Bc_(n+1,:);
    bsn = Bs_(n+1,:);

    % add to output fields
    w = w + wn.*sin(mn*zplot);
    u = u + un.*cos(mn*zplot);
    v = v + vn.*cos(mn*zplot);

    b = b + bcn.*cos(mn*zplot) + bsn.*sin(mn*zplot);
end

w = w + w0.*(1-zplot/H);
u = u + u0;
v = v + V_(1,:);

b = b + Bc_(1,:) + Bs_(1,:).*(1-zplot/H);

% wave energy density
E = 1/2*(u.^2 + v.^2 + w.^2 + b.^2/N^2);

clear wn wn_Int un vn bcn bsn W_ W_Int_ U_ V_ Bs_ Bc_ M R Wk_ W_Intk_ W_Int2k_ ...
    Uk_ Vk_ Bck_ Bsk_ Mv Mbc Mbs Rbs

%%
figure;

xCutoff = x(end);
plotInds = 1:length(x(x <= xCutoff));

field_to_plot = E;

PLOT = pcolor(xplot(:,plotInds), zplot(:,plotInds), field_to_plot(:,plotInds));
set(PLOT, 'EdgeColor', 'none');
title("\(w^\prime\) for \(\Lambda = \) " + Lam + ", \(n_{\mathrm{max}}\) = " + nMax);
colorbar;
xlabel('\(x\)');
ylabel('\(z\)');
maxclim = max(abs(clim));
clim([-maxclim,maxclim]);
hold on
yline((abs(f)/kstar - U0)/Lam,'--','LineWidth',2,'Color',[0,158,115]/255);   % critical levels
yline((-abs(f)/kstar - U0)/Lam,'--','LineWidth',2,'Color',[0,158,115]/255);
yline((N/kstar - U0)/Lam,'-.k');        % turning level
yline(-U0/Lam,'--g');
hold off

xrange = xlim;
pbaspect([xrange(2)-xrange(1),H,1]);
