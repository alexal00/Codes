% initialize variables (initial guess)
T_D0 = W; % N
H_D0 = W/20; % N
R_f0 = W/20; % N
M_f0 = 0.; % Nm
P_c0 = 1000.; % N
a_10 = 0.; % radian
theta_00 = 8/180*pi; % radian
B_10 = 0.; % radian
gamma0 = 0.; % radian
tau0 = 0.; 
V_infty0 = 50;
mu0 = V_infty0/v_tip;
lambda0 = 0.02;
alpha_D0 = 0.; % radian
u0 = lambda0*v_tip; % m/s
v10 = v_tip*sqrt(mu0^2 + lambda0^2); % m/s