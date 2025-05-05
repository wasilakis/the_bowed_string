% bowed-string simulation using an elasto-plastic friction model
% passivity is guaranteed when s1 is velocity-dependent

% Date: 05/05/2025
% Autors: Ewa Matusiak & Vasileios Chatziioannou (mdw)
% license: GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

%%%%%%%%%%%%%%%%%%%% Related publication %%%%%%%%%%%%%%%%%%%%
% Matusiak, Chatziioannou & van Walstijn,
% Numerical modelling of elasto-plastic friction in bow-string interaction with guaranteed passivity
% Front. Signal Process. (5):2025, doi: 10.3389/frsip.2025.1525044

clear;

str_par = [0.7; 5e-4; 98; 149.7415; 1.37e10; 1.53714; 0.0087];
% str_par = [L, r, f_0, T, E, sig0, sig1]  ---> string parameters
%           L = length in [m]
%           r = radius in [m]
%           f_0 = fundamental frequency in [Hz]
%           T = tension in [N]
%           E = Young's modulus in [Pa]
%           sig0 = freq. independent damping in [1/s]
%           sig1 = freq. dependent damping in [m^2/s]

str_tor_par = [3.0312e-4; 4.2e-10; 29];
% str_tor_par = [cT, theta, Q]  ---> string parameters for torsional motion
%               cT = torsional wave speed in [m/s]
%               theta = polar moment of interia
%               Q = factor for torsional waves, used for damping

bow_par = [4; 0.01];
% bow_par = [M, bow_width] ---> bow parameters
%           M = number of points to consider under the bow
%           bow_width = width of the bow in [m]

BH_par = [4.5e-3; 4.8297e4; 5.7674];
% BH_par = [bow_hair mass, spr_const, damp_const] ---> lumped bow-hair parameters
%          bow_hair mass = mass in [kg]
%          spr_const = spring constant in [N/m]
%          damp_const = damping constant in [kg/s]

ss_LuGre = [3.1860e5; 0.0027; 0; 0.2280; 0.5071; 1.0207; 2];
% ss_LuGre = [s0, s1, s2, vs, muC, muS, Sexp] ---> friction parameters
%            s0 = stiffness of the bristles in [N/m] 
%            s1 = damping of the bristles in [kg/s] 
%            s2 = viscous friction coefficient 
%            vs = Stribeck velocity
%            muC = Coulumb friction coefficient
%            muS = Stribeck friction coefficient
%            Sexp = Stribeck exponent

% slFunction = 0;             % constant s1 (no passivity proof)
s1Function = 1;             % velocity dependent s1 (guaranteed passive)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur = 0.5;              % simulation duration
fs = 44100;             % sampling rate
dt = 1/fs;              % time step

v_b = 0.3439;           % bow velocity
fN = 2.3433;            % bow force
beta = 0.0786;          % bow position (ratio over string length)
acc = 0.8722;           % bow acceleration

time = dt:dt:dur;       % time vector

TS = ceil(fs*dur);      % number of discrete time steps

 if acc > 0
    acc_dur = abs(v_b)/acc; % if smaller than dur: constant velocity
    TS_a = ceil(acc_dur*fs);
    v_b_vec = [0:v_b/(TS_a-1):v_b v_b*ones(1,TS-TS_a)];
elseif acc == 0
    v_b_vec = v_b*ones(1,TS);
else
    disp('acc can not be negative')
end


maxiter = 100;          % max number of iterations in Newton method
thrsh = 1e-15;          % threshold

order = 3;              % interpolation order

M = bow_par(1);         % desired number of points under the bow
bow_width = bow_par(2); % width of the bow in [m]

%% reading model parameters
L = str_par(1);
r = str_par(2);
f_0 = str_par(3);
T = str_par(4);
E = str_par(5);
sig0 = str_par(6);
sig1 = str_par(7);

s0 = ss_LuGre(1);
s1_bar = ss_LuGre(2);
s2 = ss_LuGre(3);
vS  = ss_LuGre(4);

muC = ss_LuGre(5);
muS = ss_LuGre(6);
Sexp = ss_LuGre(7);

%% calculated parameters

fC = muC*fN;
fS = muS*fN;

z_ba = 0.7*fC/s0;       % bristle breaking point

x_b = L*beta;           % location of the bow

% assume that bowLocS_l > h and bowLocS_r < L-h
bowLocS_l = x_b - bow_width/2;    % left end location of the bow on the string
bowLocS_r = x_b + bow_width/2;    % right end location of the bow on the string

if bow_width == 0
    bowLocS = x_b;          % chosen points under the bow to compute relative velocity
    M = length(bowLocS);    % number of points chosen under the bow for relative velocity
else
    if M == 1
        bowLocS = x_b;
    else
        h2 = bow_width/M;
        bowLocS = bowLocS_l:h2:bowLocS_r;   % chosen points under the bow to compute relative velocity
        M = length(bowLocS);                % number of points chosen under the bow for relative velocity
    end
end

%% create matrices for the string motion

% derived parameters
A = pi*r^2;           % cross sectional area (m^2)
I = pi*r^4/4;         % area moment of inertia (m^4)
c = f_0*2*L;
rho = T/(c^2*A);

dens = rho*A;         % [kg/m]
bend_stiff = E*I;     % [N/m^2]
kap = sqrt(bend_stiff/dens);   % stiffness coefficent (m^2/s)

hh = (c*dt)^2 + 4*dt*sig1;
h = sqrt((hh+sqrt(hh^2+16*(dt*kap)^2))/2);
N = floor(L/h);       % number of discrete spatial steps
h = L/N;              % length step

al = c*dt/h;
be = kap*dt/h^2;
ga = 2*sig1*dt/h^2;

sig = 1 + sig0*dt;

O1 = 2/dt * sig;

Dxx = sparse(toeplitz([-2 1 zeros(1,N-3)]));
Dxxxx = Dxx*Dxx;

Id = sparse(eye(N-1));

Am = (2*Id +(ga+al^2)*Dxx - be^2 * Dxxxx);
Bm = ((sig0*dt-1)*Id - ga*Dxx);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Iu,Ju] = interpolation(bowLocS,h,N-1,order);
norJ_mat = 1/M * Iu*Ju;

% for energy analysis
H_r = zeros(1,TS); % energy transverse
H_w = zeros(1,TS); % energy torsional
H_h = zeros(1,TS); % total energy bow
H_b = zeros(1,TS); % total bristle energy

H_t = zeros(1,TS);  % total energy
dH_t = zeros(1,TS);
h_e = zeros(1,TS);  % energy conservation error
Loss = zeros(1,TS);
%%%%%%%%%%%%%%%%%%%%%%%
Q_r = zeros(1,TS); % damping loss of transverse motion
Q_w = zeros(1,TS); % damping loss of torsional motion
Q_h = zeros(1,TS); % loss in bow hair
Q_b = zeros(1,TS); % loss in bristles

P = zeros(1,TS);
qS = 0;

%% initialization

% simply supported boundary conditions
u_p = zeros(N-1,1);   % time zero
u = u_p;     % time one

v = zeros(M,TS); % relative velocity

F_fr = zeros(M,TS);
F_bridge = zeros(1,TS);
F_nut = zeros(1,TS);

% vectors for testing
Niter = zeros(1,TS); % Newton iterations

% the main update loop
v(:,1) = -v_b_vec(1);
vr = v(:,1);

%%%% friction update %%%
z = zeros(M,1);
zmh = zeros(M,1);
z_av = 1/2 * (z + zmh);

% compute s1
if s1Function == 1
    pow = 1;
    pp = 2*pow;
    epsilon = fC/s1_bar;
    denom = (vr.^pp + epsilon^pp).^(1/pp);
    s1 = fC./denom;
end

if s1Function == 0
    s1 = s1_bar;
    d_s1 = 0;
end

% create matrices for the torsional string motion
K_T = str_tor_par(1);
P_T = str_tor_par(2);
Qfac = str_tor_par(3);
cT = sqrt(K_T/P_T);

%%%% derived parameters
if Qfac == 0
    sig2 = 0;
else
    sig2 = 1/(2*Qfac);        % damping coefficient
end

sigQ = sig2;
sigT = 1+sigQ*dt;

O2 = 2/dt*sigT;

hT = cT*dt;
NT = floor(L/hT);        % number of discrete spatial steps for torsional motion
hT = L/NT;               % length step torsional

muu = cT^2 * dt^2/hT^2;
zeta = r*dt^2/P_T;

DxxT = sparse(toeplitz([-2 1 zeros(1,NT-3)]));         % acting on sol(2:N)

IdT = sparse(eye(NT-1));

Cm = (2*IdT + muu*DxxT);

% fixed boundary conditions
w_p = zeros(NT-1,1);         % time zero
w = w_p;           % time one

% bow hair
dens_h = BH_par(1);
spr_const = BH_par(2);
damp_const = BH_par(3);

% lumped -> distributed parameters
if bow_width > 0
    dens_h = dens_h/bow_width;
    spr_const = spr_const/bow_width;
    damp_const = damp_const/bow_width;
end

c_h1 = dens_h /(dt^2) + spr_const/4 + damp_const/(2*dt);
c_h2 = dens_h /(dt^2) + spr_const/4 - damp_const/(2*dt);

eta = zeros(M,1);
eta_p = eta;

O3 = 2/dt*dens_h + dt/2 * spr_const + damp_const;

[IuT,JuT] = interpolation(bowLocS,hT,NT-1,order);
norJT_mat = 1/M * IuT*JuT;

IJ_mat = 1/dens*1/O1*norJ_mat + h/hT * r^2 * 1/P_T * 1/O2 *norJT_mat + 1/M * 1/O3;

for n = 2:TS

    vB = v_b_vec(n);

    % auxiliary variables
    d_u = (u - u_p)/dt;
    d_w = (w - w_p)/dt;
    d_eta = (eta - eta_p)/dt;

    su = (c^2/O1 * 1/h^2 * Dxx - kap^2/O1 * 1/h^4 * Dxxxx)*u + 1/O1 *(2/dt * Id + 2/h^2 * sig1*Dxx)*d_u;
    sun = Iu*su;

    sw = r*(cT^2/O2 * 1/hT^2 *DxxT*w + 2/dt * 1/O2 *d_w);
    swn = h/hT * IuT*sw;

    dn = ((dt/2 * spr_const + 2/dt * dens_h)*d_eta - spr_const*eta)/O3;

    sn = sun - swn + dn - vB;

    err = 1;
    i = 0;
    Newton = true;

    while Newton
        % compute s1 and derivative of s1
        if s1Function == 1
            denom = (vr.^pp + epsilon^pp).^(1/pp);
            s1 = fC./denom;
            d_s1 = - fC *vr.^(pp-1).*(1./denom).^(pp+1);
        end

        dz = 2/dt * (z_av - zmh);

        % compute z_ss & derivative wrt vr
        z_ss = sign(vr).*(fC + (fS - fC)*exp(-(vr/vS).^Sexp) + s2*abs(vr))./s0;
        dz_ss = ((-Sexp*abs(vr).^(Sexp-1)).*(sign(vr).^Sexp)./(vS^Sexp*s0)).*((fS-fC)*exp(-(vr/vS).^Sexp)) + sign(vr)*s2./s0;

        ind = vr == 0;
        z_ss(ind) = fS./s0;
        dz_ss(ind) = 2*fS./s0;

        % compute derivatives of theta & alpha

        alpha = 1/2*(1+sign(z_av).*sin(pi*(z_av - 1/2*sign(z_av).*(abs(z_ss)+z_ba))./(abs(z_ss)-z_ba)));

        dz_ss_vAbs = sign(z_ss).*dz_ss;

        theta = pi*(abs(z_av)-(abs(z_ss)+z_ba)/2)./(abs(z_ss)-z_ba);

        d_alpha_z = sign(z_av).*(pi/2*cos(sign(z_av).*theta)./(abs(z_ss)-z_ba));
        d_alpha_v = dz_ss_vAbs.*((z_ba-abs(z_av))./(abs(z_ss)-z_ba).^2*pi/2.*cos(sign(z_av).*theta));

        ind1 = (sign(vr) == sign(z_av) & abs(z_av) <= z_ba) | sign(vr) ~= sign(z_av);
        alpha(ind1) = 0;
        d_alpha_z(ind1) = 0;
        d_alpha_v(ind1) = 0;

        ind2 = sign(vr) == sign(z_av) & abs(z_av) >= abs(z_ss);
        alpha(ind2) = 1;
        d_alpha_z(ind2) = 0;
        d_alpha_v(ind2) = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        gn = vr.*(1 - alpha.*z_av./z_ss);

        fvz = s0*z_av + s1.*gn;

        N1 = vr + IJ_mat*fvz - sn;
        N2 = gn - dz;

        dgn_z = -vr./z_ss.*(d_alpha_z.*z_av + alpha);
        dgn_v = 1 - z_av.*((alpha + d_alpha_v.*vr).*z_ss - dz_ss.*alpha.*vr)./z_ss.^2;

        % Jacobian
        dN1v = eye(M) + IJ_mat*diag(s1.*dgn_v + d_s1.*gn);
        dN1z = IJ_mat*diag(s0 + dgn_z.*s1);
        dN2v = diag(dgn_v);
        dN2z = diag(dgn_z - 2/dt);

        dN = [dN1v dN1z; dN2v dN2z];

        % increase Jacobian when convergence is slow (optional)
        if i > 50
            dN = 1.1*dN;
        end

        % Newton Raphson iteration step
        vzPrev = [vr; z_av];
        vzNew = vzPrev - dN\[N1;N2];

        err = norm(vzNew-vzPrev);

        vr = vzNew(1:M,1);
        z_av = vzNew(M+1:2*M,1);

        i = i+1;
        Newton = (err > thrsh && i < maxiter);

    end
    %%%%%%%% end of Newton's method

    if s1Function == 1
        denom = (vr.^pp + epsilon^pp).^(1/pp);
        s1 = fC./denom;
    end

    zph = 2*z_av - zmh;
    gn = (zph - zmh)/dt;

    F_fr(:,n) = s0*z_av + s1.*gn;

    u_n = 1/sig*(Am*u + Bm*u_p) - dt^2*Ju*F_fr(:,n)/(M*sig*dens);
    w_n = 1/sigT *(Cm*w - (1 - sigQ*dt)*w_p + zeta*JuT*F_fr(:,n)/M);
    eta_n = ((2*dens_h/(dt^2) - spr_const/2)*eta - c_h2*eta_p - F_fr(:,n)/M)/c_h1;

    % Energy
    vrrr = Iu*((u_n - u_p)/(2*dt)) - r * h/hT * IuT*((w_n - w_p)/(2*dt)) + (eta_n - eta_p)/(2*dt) - vB;

    H_r(n) = dens*h/2*sum((1/dt*(u - u_p)).^2) + dens*(c^2/2 * 1/h *sum(([u;0] - [0;u]).*([u_p;0] - [0;u_p])) + kap^2/(2*h^3)*sum((Dxx*u).*(Dxx*u_p)));
    H_w(n) = h/hT * (P_T* hT/2*sum((1/dt)^2*(w - w_p).^2) + K_T/(2*hT)*sum(([w;0] - [0;w]).*([w_p;0] - [0;w_p])));
    H_h(n) = dens_h/2*sum((1/dt*(eta - eta_p)).^2) + spr_const/2 * sum(((eta + eta_p)/2).^2);
    H_b(n) = 1/M * s0/2 * sum(zmh.^2);

    Q_r(n) = dens*(h*sig0*sum((u_n - u_p).^2)/(2*dt^2)- sig1/(h*dt^2)*sum((u_n - u_p).*(Dxx*(u - u_p))));
    Q_w(n) = h/hT*P_T*hT*sigQ*sum((w_n - w_p).^2)/(2*dt^2);
    Q_h(n) = damp_const*sum(((eta_n - eta_p)/(2*dt)).^2);
    Q_b(n) = 1/M * sum(vrrr.*F_fr(:,n)) - 1/M *sum(s0*z_av.*gn);

    P(n) = 1/M * sum(vB.*F_fr(:,n));

    H_t(n) = H_r(n) + + H_w(n) + + H_h(n) + H_b(n);
    Loss(n) = Q_r(n) + Q_w(n) + Q_h(n) + Q_b(n) + P(n);
    h_e(n) = (H_t(n) + qS - H_t(1));

    qS = qS+dt*Loss(n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the variables
    zmh = zph;

    u_p = u;
    u = u_n;

    w_p = w;
    w = w_n;

    eta_p = eta;
    eta = eta_n;

    % store vectors
    v(:,n) = vr;
    Niter(n) = i;

    F_bridge(n) = dens*(c^2/h * u_n(1) - kap^2/h^3 * (u_n(2) - 2*u_n(1)));
    F_nut(n) = dens*(c^2/h * u_n(N-1) - kap^2/h^3 * (u_n(N-2) - 2*u_n(N-1)));
end

%% plotting

figure();
ax(1) = subplot(211);
plot(time,F_bridge)
xlabel('time [s]')
ylabel('bridge force [N]')
ax(2) = subplot(212);
plot(time,h_e,'.')
xlabel('time [s]')
ylabel('energy error')

%% end of code (interpolation function below)

function [Iu,Ju] = interpolation(x,h,NS,order)

% Lagrange Interpolation Operator at point x = [x_1,...,x_M]
% h = grid spacing
% NS = number of grid points (N+1);
% order of the interpolation = 0,1 or 3;
N = NS-1;
M = length(x);

Iu = zeros(M,NS);        % interpolation operator

if order == 0
    for i=1:M
        xcu = floor(x(i)/h);
        %alphaU = x(i)/h - xcu;
        xcu = xcu + 1; %MATLAB indexind starting from 1
        Iu(i,xcu) = 1;
        if xcu == N+2
            Iu(i,xcu-1) = 1;
        end
    end
    
elseif order == 1
    
    for i=1:M
        xcu = floor(x(i)/h);
        alphaU = x(i)/h - xcu;
        xcu = xcu + 1; %MATLAB indexind starting from 1
            
        if xcu == N+1
           Iu(i,xcu) = 1;
        elseif xcu == N+2
            Iu(i,xcu-1) = 1;
        else
           Iu(i,xcu) = 1 - alphaU;
           Iu(i,xcu+1) = alphaU;
        end
    end

elseif order == 3
     
    for i=1:M
        xcu = floor(x(i)/h);
        alphaU = x(i)/h - xcu;
        xcu = xcu + 1;       % MATLAB indexind starting from 1
 
        if xcu == 1
           Iu(i,xcu) = (alphaU - 1)*(alphaU + 1)*(alphaU - 2)/2;
           Iu(i,xcu+1) = -alphaU*(alphaU + 1)*(alphaU - 2)/2 + alphaU*(alphaU - 1)*(alphaU - 2)/6;
           Iu(i,xcu+2) = alphaU*(alphaU + 1)*(alphaU - 1)/6; 
        elseif xcu == N
           Iu(i,xcu-1) = -alphaU*(alphaU - 1)*(alphaU - 2)/6;
           Iu(i,xcu) = (alphaU - 1)*(alphaU + 1)*(alphaU - 2)/2 - alphaU*(alphaU + 1)*(alphaU - 1)/6;
           Iu(i,xcu+1) = -alphaU*(alphaU + 1)*(alphaU - 2)/2;
        elseif xcu == N+1
           Iu(i,xcu) = (alphaU - 1)*(alphaU + 1)*(alphaU - 2)/2;
        elseif xcu == N+2
            Iu(i,xcu-1) = 1;
        else
           Iu(i,xcu-1) = -alphaU*(alphaU - 1)*(alphaU - 2)/6;
           Iu(i,xcu) = (alphaU - 1)*(alphaU + 1)*(alphaU - 2)/2;
           Iu(i,xcu+1) = -alphaU*(alphaU + 1)*(alphaU - 2)/2;
           Iu(i,xcu+2) = alphaU*(alphaU + 1)*(alphaU - 1)/6;
        end
    end
end

Ju = 1/h * Iu';      % spreading operator

end

