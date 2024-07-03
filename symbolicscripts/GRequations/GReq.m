close all; clear variables; clc
%% Symbolic development of the GR equations
% The GR equations of motion will be derived from the simplified model
% described in Meirovitch, Coleman and Hammond.

% Airframe parameters
% X-direction parameters
syms Mx Cx Kx real
% Y-direction parameters
syms My Cy Ky real
% Blade parameters
syms Se Nb M Om Je Cd Ced Ke e real

% Equations of motion

% Single blade equation of motion in the rotation reference
blades = 5;
% Mass matrix
MRb = Je*eye(blades);

% Damping matrix
CRb = sym(Cd*eye(blades));
C = [Ced,0,Cd,0,Ced];

% Populate the sub-diagonals in a circular manner
for i = 1:blades
    % Calculate the indices wrapping around
    idx1 = mod(i-3, blades) + 1; % i-2 with wrapping
    idx2 = mod(i+1, blades) + 1;  % i+2 with wrapping
    
    % Populate the sub-diagonals
    CRb(i, idx1) = CRb(i, idx1)-Ced; % Populate a(i,i-2)
    CRb(i, idx2) = CRb(i, idx2)-Ced; % Populate a(i,i+2)
end

% CRb = subs(CRb,Ced,Cd);

% Stiffness matrix
kd = Ke+Om^2*e*Se;
KRb = kd*eye(blades);
% Lead-lag rotating coordinates
qR = sym('chi',[blades 1]);

% Mutiblade coordinates
T = sym(zeros(blades));
syms chi_0
if mod(blades,2)==0
    chi_S = sym(['chi_' num2str(blades/2)]);
    chi_c = sym('chi_c',[blades/2-1 1]);
    chi_s = sym('chi_s',[blades/2-1 1]);
else
    chi_c = sym('chi_c',[(blades-1)/2 1]);
    chi_s = sym('chi_s',[(blades-1)/2 1]);
end

qNR = sym(zeros(blades,1));
qNR(1) = qNR(1)+chi_0;
aux = 2;
for ii=1:length(chi_c)
    qNR(aux) = chi_c(ii);
    qNR(aux+1) = chi_s(ii);
    aux = aux + 2;
end
if exist("chi_S","var")
    qNR(end) = chi_S(1);
end
% syms deltapsi real
deltapsi = 2*sym(pi)/blades;

psi = sym(zeros(blades,1));
T = sym(zeros(blades));
Sb_R = sym(zeros(blades,2));
for ii = 1:blades
    psi(ii) = sym('psi1')+(ii-1)*deltapsi;
    Sb_R(ii,[1 2])=Se*[-sin(psi(ii)) cos(psi(ii))];
    T(:,1) = 1;
    aux = 2;
    for jj=1:length(chi_c)
    T(ii,aux) = cos((jj)*psi(ii));
    T(ii,aux+1) = sin((jj)*psi(ii));
    aux = aux+2;
    end
    if exist("chi_S","var")
    T(ii,end) = (-1)^ii;
    end
end

% Derivation of matrices wrt time
syms t

% Inertial coupling of matrices
Sb_R = subs(Sb_R,psi(1),Om*t);
Sb_Rdot = diff(Sb_R,t); Sb_Rdot = subs(Sb_Rdot,t,0);
Sb_Rddot = diff(Sb_R,t,2); Sb_Rddot = subs(Sb_Rddot,t,0);
Sb_R = subs(Sb_R,t,0);
% Rotation matrix
T = subs(T,psi(1),Om*t);
Tdot = diff(T,t); Tdot = simplify(subs(Tdot,t,0));
Tddot = diff(T,t,2); Tddot = simplify(subs(Tddot,t,0));
T = simplify(subs(T,t,0));

%% Equations of motion in the NR frame
for ii=1:length(chi_c)
    Epsilon(ii) = 2*(1-cos(2*ii*deltapsi));
end
% Blade equations
% Mass matrix
MbNRb = T'*MRb*T; MbNRb = simplify(MbNRb);
MbNRb = subs(MbNRb,blades,Nb);
if mod(blades,2)==0
   MbNRb = subs(MbNRb,blades/2,Nb/2);
end
MbNRa = T'*Sb_R; MbNRa = simplify(MbNRa);
% MbNRa = subs(MbNRa,blades,Nb);
% MbNRa = subs(MbNRa,-blades,-Nb);
if mod(blades,2)==0
   MbNRa = subs(MbNRa,blades/2,Nb/2);
   MbNRa = subs(MbNRa,-blades/2,-Nb/2);
end
% Damping matrix
CbNRb = T'*(CRb*T+2*MRb*Tdot); CbNRb = simplify(CbNRb);
% CbNRb = subs(CbNRb,blades,Nb)
% if mod(blades,2)==0
%    CbNRb = subs(CbNRb,blades/2,Nb/2);
%    CbNRb = subs(CbNRb,-blades/2,-Nb/2);
% end

% Stiffness matrix
KbNRb = T'*(KRb*T+CRb*Tdot+MRb*Tddot); KbNRb = simplify(KbNRb,"Steps",50);
%KbNRb = KbNRb/5; %KbNRb = KbNRb*5;

% Airframe equations
Ma = [Mx+Nb*M,          0;...
            0,    My+Nb*M];
Ca = [     Cx,          0;...
            0,         Cy];
Ka = [     Kx,          0;...
            0,          Ky];
% Mass matrix
MaNRa = Ma;
CaNRa = Ca;
KaNRa = Ka;
MaNRb = Sb_R'*T; MaNRb = simplify(MaNRb);

dof_col = 1;
dof_cou = [2 3];

if length(chi_c)>1
    aux = 4;
    aux2 = 1;
    for ii = 1:length(chi_c)-1
    dof_cyc(aux2) = aux;
    dof_cyc(aux2+1) = aux+1;
    aux=aux+2;
    aux2 = aux2+2;
    end
end

% dof
COLM = constmat(MbNRb,MbNRa,zeros(2),dof_col);
COLC = constmat(CbNRb,zeros(6,2),zeros(2),dof_col);
COLK = constmat(KbNRb,zeros(6,2),zeros(2),dof_col);

GRM = constmat(MbNRb,MbNRa,MaNRa,dof_cou);
GRC = constmat(CbNRb,zeros(6,2),CaNRa,dof_cou);
GRK = constmat(KbNRb,zeros(6,2),KaNRa,dof_cou);
% COLC
if exist("chi_S")
    dof_sci = blades;
    SCIM = constmat(MbNRb,MbNRa,zeros(2),dof_sci);
    SCIC = constmat(CbNRb,zeros(6,2),zeros(2),dof_sci);
    SCIK = constmat(KbNRb,zeros(6,2),zeros(2),dof_sci);
end


