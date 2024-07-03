%% Definiton of the basic helicopter parameters
% Data for verification:
% helicopter basic parameters
R = 6.06; data.R = R;
Omega = 337.8*2*pi/60; data.Omega =Omega;
b = 5; data.b = b;
e = 0.04*R; data.e = e;
Mb = 70.; data.Mb = Mb;
% uniform blade assumption
Sb = Mb*(R - e)/2;
data.Sb = Sb;
Ib = Mb*(R - e)^2/3; 
data.Ib = Ib;
K_beta = 0.; data.Kb = K_beta;

data.y = blade(R,e,100);
data.nu_b = sqrt(1+e*Sb/Ib+K_beta/(Ib*Omega^2));

% global parameters (used later to compute function and Jacobian matrix)
% global h xcg l Kh A v_tip sigma cla cd0 f Sbar lbar cmf clt St W rho
h = 2.; data.h =h; % m
xcg = 0.; data.xcg = xcg; % m
l = 12.; data.l = l;% m
Kh = b/2*(e*Sb*Omega^2 + K_beta); data.Kh = Kh;% N m/radian
A = pi*R^2; data.A = A;% m^2
v_tip = R*Omega; data.v_tip=v_tip;% m/s
sigma = 0.08; data.sigma = sigma;
c = pi*R*sigma/b; data.c = c;
B = 1-c/(2*R); data.B = B;       % Tip loss factor
data.thetaw = 0;                % blade linear twist
% C_Lalpha = 5.73; % typical coefficient (91% of 2*pi)
% C_Lalpha = 6.0447; % NACA 23012 at Mach 0-0.2 used in MBDyn
cla = 7.0760; data.cla=cla;% NACA 23012 at Mach 0.5 used in MBDyn
cd0 = 0.008; data.cd0=cd0;
f = .4; data.f =f;% m^2
Sbar = 2.; data.Sbar=Sbar;% m^2
lbar = 1.; data.lbar = lbar;% m
cmf = 0.02; data.cmf=cmf;
clt = 5.73; data.clt=clt;

St = 1.2; data.St=St;% m^2
W = 45200.; data.W=W;% N
rho = 1.225; data.rho=rho;% kg/m^3
a = 340.3; data.a = a; % m/s

