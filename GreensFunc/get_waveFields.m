function [w,u,v,b,phi,...
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
                                w0_Int,w0_Int2,w0_Int3,zplot,nMax,f,N,U0,H,Lam)

% initialize fields
w = 0*zplot;
w_x = 0*zplot;
w_x2 = 0*zplot;
w_x3 = 0*zplot;
w_Int = 0*zplot;
w_Int2 = 0*zplot;

w_z = 0*zplot;
w_xz = 0*zplot;
w_x2z = 0*zplot;
w_x3z = 0*zplot;
w_Int_z = 0*zplot;
w_Int2_z = 0*zplot;
w_Int3_z = 0*zplot;

w_z2 = 0*zplot;
w_xz2 = 0*zplot;
w_x2z2 = 0*zplot;
w_Int_z2 = 0*zplot;
w_Int2_z2 = 0*zplot;
w_Int3_z2 = 0*zplot;

w_z3 = 0*zplot;
w_xz3 = 0*zplot;
w_Int_z3 = 0*zplot;
w_Int2_z3 = 0*zplot;
w_Int3_z3 = 0*zplot;


% modal parts -------------------------------------------------------------
for n = 1:nMax
    mn = m_(n);

    wn = W_(n,:);
    wn_x = W_(nMax+n,:);
    wn_x2 = W_(2*nMax+n,:);
    wn_x3 = W_(3*nMax+n,:);
    wn_Int = W_Int3_(2*nMax+n,:);
    wn_Int2 = W_Int3_(nMax+n,:);
    wn_Int3 = W_Int3_(n,:);

    % add to output field
    % w and its x-derivatives and integrals
    w = w + wn.*sin(mn*zplot);
    w_x = w_x + wn_x.*sin(mn*zplot);
    w_x2 = w_x2 + wn_x2.*sin(mn*zplot);
    w_x3 = w_x3 + wn_x3.*sin(mn*zplot);
    w_Int = w_Int + wn_Int.*sin(mn*zplot);
    w_Int2 = w_Int2 + wn_Int2.*sin(mn*zplot);

    % w_z and its x-derivatives and integrals
    w_z = w_z + mn*wn.*cos(mn*zplot);
    w_xz = w_xz + mn*wn_x.*cos(mn*zplot);
    w_x2z = w_x2z + mn*wn_x2.*cos(mn*zplot);
    w_x3z = w_x3z + mn*wn_x3.*cos(mn*zplot);
    w_Int_z = w_Int_z + mn*wn_Int.*cos(mn*zplot);
    w_Int2_z = w_Int2_z + mn*wn_Int2.*cos(mn*zplot);
    w_Int3_z = w_Int3_z + mn*wn_Int3.*cos(mn*zplot);

    % w_z2 and its x-derivatives and integrals
    w_z2 = w_z2 + -mn^2*wn.*sin(mn*zplot);
    w_xz2 = w_xz2 + -mn^2*wn_x.*sin(mn*zplot);
    w_x2z2 = w_x2z2 + -mn^2*wn_x2.*sin(mn*zplot);
    w_Int_z2 = w_Int_z2 + -mn^2*wn_Int.*sin(mn*zplot);
    w_Int2_z2 = w_Int2_z2 + -mn^2*wn_Int2.*sin(mn*zplot);
    w_Int3_z2 = w_Int3_z2 + -mn^2*wn_Int3.*sin(mn*zplot);

    % integrals of w_z3
    w_z3 = w_z3 + -mn^3*wn.*cos(mn*zplot);
    w_xz3 = w_xz3 + -mn^3*wn_x.*cos(mn*zplot);
    w_Int_z3 = w_Int_z3 + -mn^3*wn_Int.*cos(mn*zplot);
    w_Int2_z3 = w_Int2_z3 + -mn^3*wn_Int2.*cos(mn*zplot);
    w_Int3_z3 = w_Int3_z3 + -mn^3*wn_Int3.*cos(mn*zplot);

end

% not modal part ---------------------------------------------------------
a1 = -w0;
% get derivatives and integrals of a1
a1_x = -w0_x;
a1_x2 = -w0_x2;
a1_x3 = -w0_x3;

a1_Int = -w0_Int;
a1_Int2 = -w0_Int2;
a1_Int3 = -w0_Int3;


w = w + ( w0 + a1.*(zplot/H) );
w_x = w_x + ( w0_x + a1_x.*(zplot/H) );
w_x2 = w_x2 + ( w0_x2 + a1_x2.*(zplot/H) );
w_x3 = w_x3 + ( w0_x3 + a1_x3.*(zplot/H) );
w_Int = w_Int + ( w0_Int + a1_Int.*(zplot/H) );
w_Int2 = w_Int2 + ( w0_Int2 + a1_Int2.*(zplot/H) );

w_z = w_z + ( (1/H)*a1 );
w_xz = w_xz + ( (1/H)*a1_x );
w_x2z = w_x2z + ( (1/H)*a1_x2 );
w_x3z = w_x3z + ( (1/H)*a1_x3 );
w_Int_z = w_Int_z + ( (1/H)*a1_Int );
w_Int2_z = w_Int2_z + ( (1/H)*a1_Int2 );
w_Int3_z = w_Int3_z + ( (1/H)*a1_Int3 );

% higher derivatives of w
w_x4 = 2*f^2*Lam*w_z./(U0 + Lam*zplot).^3 - N^2*w_x2./(U0 + Lam*zplot).^2 ...
        - f^2*w_z2./(U0 + Lam*zplot).^2 - w_x2z2./(U0 + Lam*zplot).^3;


% fields that aren't w ----------------------------------------------------
u = -w_Int_z;
v = f*w_Int2_z./(U0 + Lam*zplot);
b = -N^2*w_Int./(U0 + Lam*zplot) + f^2*Lam*w_Int3_z./(U0 + Lam*zplot).^2;
phi = (U0 + Lam*zplot).*w_Int_z + f^2*w_Int3_z./(U0 + Lam*zplot) - Lam*w_Int;

xi = Lam*w_Int2./(U0 + Lam*zplot).^2 - w_Int2_z./(U0 + Lam*zplot);
eta = f*w_Int3_z./(U0 + Lam*zplot).^2;

zeta = w_Int./(U0 + Lam*zplot);
zeta_z = (w_Int_z.*(U0 + Lam*zplot) - Lam*w_Int)./(U0 + Lam*zplot).^2;

u_z = -w_Int_z2;
v_z = f*( w_Int2_z2.*(U0 + Lam*zplot) - Lam*w_Int2_z )./(U0 + Lam*zplot).^2;
b_z = -N^2*( w_Int_z.*(U0 + Lam*zplot) - Lam*w_Int )./(U0 + Lam*zplot).^2 ...
        + f^2*Lam* ( w_Int3_z2.*(U0 + Lam*zplot).^2 - 2*Lam*(U0 + Lam*zplot).*w_Int3_z )./(U0 + Lam*zplot).^4;
phi_z = Lam*w_Int_z + (U0 + Lam*zplot).*w_Int_z2 ...
        + f^2*( w_Int3_z2.*(U0 + Lam*zplot) - Lam*w_Int3_z )./(U0 + Lam*zplot).^2 ...
        - Lam*w_Int_z;

u_z2 = -w_Int_z3;
v_z2 = f*w_Int2_z3./(U0 + Lam*zplot) - 2*f*Lam*w_Int2_z2./(U0+Lam*zplot).^2 ...
        + 2*f*Lam^2*w_Int2_z./(U0+Lam*zplot).^3;
b_z2 = -N^2*w_Int_z2./(U0 + Lam*zplot) + 2*N^2*Lam*w_Int_z./(U0 + Lam*zplot).^2 ...
        - 2*N^2*Lam^2*w_Int./(U0 + Lam*zplot).^3 + 2*f^2*Lam*w_Int3_z3./(U0 + Lam*zplot).^2 ...
        - 4*f^2*Lam^2*w_Int3_z2./(U0 + Lam*zplot).^3 ...
        + 6*f^2*Lam^3*w_Int3_z./(U0 + Lam*zplot).^4;


% find some fields from other fields already calculated
u_x = -w_z;
v_x = -f*u./(U0 + Lam*zplot);
b_x = (-N^2*w + f*Lam*v)./(U0 + Lam*zplot);

u_x2 = -w_xz;
v_x2 = -f*u_x./(U0 + Lam*zplot);
b_x2 = (-N^2*w_x + f*Lam*v_x)./(U0 + Lam*zplot);

u_xz = -w_z2;
v_xz = -f*u_z./(U0 + Lam*zplot) + f*Lam*u./(U0 + Lam*zplot).^2;
b_xz = f*Lam*v_z./(U0 + Lam*zplot) - N^2*w_z./(U0 + Lam*zplot) - Lam*b_x./(U0 + Lam*zplot);

u_x3 = -w_x2z;
v_x3 = -f*u_x2./(U0 + Lam*zplot);
b_x3 = (-N^2*w_x2 + f*Lam*v_x2)./(U0 + Lam*zplot);

u_x2z = -w_xz2;
v_x2z = -f*u_xz./(U0 + Lam*zplot) + f*Lam*u_x./(U0 + Lam*zplot).^2;
b_x2z = ( f*Lam*v_xz - N^2*w_xz - Lam*b_x2 )./(U0 + Lam*zplot);

u_xz2 = -w_z3;
v_xz2 = -f*u_z2./(U0 + Lam*zplot) + -2*Lam*v_xz./(U0 + Lam*zplot);
b_xz2 = ( f*Lam*v_z2 - N^2*w_z2 - 2*Lam*b_xz )./(U0 + Lam*zplot);

u_x4 = -w_x3z;
v_x4 = -f*u_x3./(U0 + Lam*zplot);
b_x4 = (-N^2*w_x3 + f*Lam*v_x3)./(U0 + Lam*zplot);

u_x3z = -w_x2z2;
v_x3z = -f*u_x2z./(U0 + Lam*zplot) + f*Lam*u_x2./(U0 + Lam*zplot).^2;
b_x3z = ( f*Lam*v_x2z - N^2*w_x2z - Lam*b_x3 )./(U0 + Lam*zplot);

u_x2z2 = -w_xz3;
v_x2z2 = -f*u_xz2./(U0 + Lam*zplot) + -2*Lam*v_x2z./(U0 + Lam*zplot);
b_x2z2 = ( f*Lam*v_xz2 - N^2*w_xz2 - 2*Lam*b_x2z )./(U0 + Lam*zplot);