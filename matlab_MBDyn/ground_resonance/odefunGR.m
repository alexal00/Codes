function dydt = odefunGR(time,y,omega,mode,blades)
% Helicopter data
data = hammond1974data;
data.n_b = blades;
data.omega = omega;
dpsi = 2*pi/blades;
if mode ==0
    fact = 1;
elseif mode==1
    fact= 1/(2*(1-cos(dpsi)));
elseif mode==2
    fact= 1/(2*(1-cos(2*dpsi)));
end
data.c_xi = fact*data.c_xi;
% syms t real
subs_var = [data.m_x data.c_x data.k_x data.m_y data.c_y data.k_y data.S_b...
            data.n_b data.m_b data.omega data.I_b data.c_xi data.k_xi data.e];
[M_R, C_R, K_R, sym_var] = symbolicGReq(mode,blades);
M_R = subs(M_R,sym_var, subs_var);
M_R = matlabFunction(M_R);
M_R = M_R(time);
C_R = subs(C_R,sym_var, subs_var);
C_R = matlabFunction(C_R);
C_R = C_R(time);
K_R = subs(K_R, sym_var, subs_var);
K_R = matlabFunction(K_R);
K_R = K_R(time);

dof = data.n_b+2;
A = [zeros(dof) eye(dof);...
     -M_R\K_R -M_R\C_R];
dydt = A*y;
end