function [data,input,x0]=probdef(ptype)
% Obtain information of the helicopter
% hel_info;
aw169;
initialcondtrim;

if strcmp(ptype,'inverse')
    input.V = V_infty0; % m/s
    input.tau = tau0; % radian
    x0 = [T_D0 H_D0 R_f0 M_f0 P_c0 a_10 gamma0 mu0 lambda0 alpha_D0 u0 v10  theta_00 B_10];
elseif strcmp(ptype,'direct')
    input.theta0 = 5.31341*pi/180; % radian
    input.B1 = 3.35661*pi/180; % radian
    alpha_H = 0.933894*pi/180;
    tau0 = alpha_H-gamma0;
    x0 = [T_D0 H_D0 R_f0 M_f0 P_c0 a_10 gamma0 mu0 lambda0 alpha_D0 u0 v10  V_infty0 tau0];
else
    Q = 0;
    input.V = V_infty0; % m/s
    input.Q = Q; % radian
    x0 = [T_D0 H_D0 R_f0 M_f0 P_c0 a_10 gamma0 mu0 lambda0 alpha_D0 u0 v10  theta_00 B_10 tau0];
end

end