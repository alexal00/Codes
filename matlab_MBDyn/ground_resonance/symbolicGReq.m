function [M_R, C_R, K_R, sym_var] = symbolicGReq(mode,blades)
%% Symbolic development of the GR equations
% The GR equations of motion
% Airframe parameters
% X-direction parameters
syms m_x c_x k_x real
% Y-direction parameters
syms m_y c_y k_y real
% Blade parameters
syms S_b n_b m_b omega I_b c_xi c_ed k_xi e real

% Equations of motion
% Single blade equation of motion in the rotation reference
% blades = 4;
airframe = 2;
% Mass matrix
MRb = I_b*eye(blades);

% Damping matrix
CRb = sym(c_xi*eye(blades));
% C = [c_ed,0,c_xi,0,c_ed];

% Populate the sub-diagonals in a circular manner
for i = 1:blades
    % Calculate the indices wrapping around
    if mode~=0
        if mode==1
        idx1 = mod(i-2, blades) + 1; % i-1 with wrapping
        idx2 = mod(i, blades) + 1;  % i+1 with wrapping
        else
        idx1 = mod(i-3, blades) + 1; % i-2 with wrapping
        idx2 = mod(i+1, blades) + 1;  % i+2 with wrapping
        end
    % Populate the sub-diagonals
    CRb(i, idx1) = CRb(i, idx1)-c_ed; % Populate a(i,i-2)
    CRb(i, idx2) = CRb(i, idx2)-c_ed; % Populate a(i,i+2)
    end
    
end

if mode == 0
CRb = subs(CRb,c_ed,0);
else
CRb = subs(CRb,c_xi,2*c_xi);
CRb = subs(CRb,c_ed,c_xi);
end

% Stiffness matrix
kd = k_xi+omega^2*e*S_b;
KRb = kd*eye(blades);
% Lead-lag rotating coordinates
% qR = sym('chi',[blades 1]);

% syms deltapsi real
deltapsi = 2*sym(pi)/blades;

psi = sym(zeros(blades,1));
% T = sym(zeros(blades));
Sb_R = sym(zeros(blades,2));
for ii = 1:blades
    psi1 = sym('psi1');
    psi(ii) = psi1+(ii)*deltapsi;
    Sb_R(ii,[1 2])=S_b*[-sin(psi(ii)) cos(psi(ii))];
end

% Derivation of matrices wrt time
syms t real

% Inertial coupling of matrices
Sb_R = subs(Sb_R,psi1,omega*t);
Sb_Rdot = diff(Sb_R,t); %Sb_Rdot = subs(Sb_Rdot,t,0);
Sb_Rddot = diff(Sb_R,t,2); %Sb_Rddot = subs(Sb_Rddot,t,0);
% Sb_R = subs(Sb_R,t,0);


% Airframe equations
Ma = [m_x+n_b*m_b,          0;...
            0,    m_y+n_b*m_b];
Ca = [     c_x,          0;...
            0,         c_y];
Ka = [     k_x,          0;...
            0,          k_y];

M_R = [MRb Sb_R;...
       Sb_R' Ma];
C_R = [CRb zeros(blades,airframe);...
       2*Sb_Rdot' Ca];
K_R = [KRb zeros(blades,airframe);...
       Sb_Rddot' Ka];

sym_var = [m_x c_x k_x m_y c_y k_y S_b n_b m_b omega I_b c_xi k_xi e];
end

