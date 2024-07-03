close all; clear all; clc
%% Declaration of variables
syms a b f d c_a c_b c_f c_d  Nb e real positive
syms xi1 xi3 deltad2 deltaf2 phi deltapsi real

%%
A = [-a;c_a;0];
Ra21 = [cos(xi1) sin(xi1) 0;...
        -sin(xi1) cos(xi1) 0;...
        0 0 1]; 
Ra10 = [cos(deltapsi) sin(deltapsi) 0;...
        -sin(deltapsi) cos(deltapsi) 0;...
        0 0 1];

A = Ra10*([0;e;0]+Ra21*A);

B = [b;c_b;0];

Rb21 = [cos(xi3) sin(xi3) 0;...
        -sin(xi3) cos(xi3) 0;...
        0 0 1]; 
Rb10 = [cos(deltapsi) -sin(deltapsi) 0;...
        sin(deltapsi) cos(deltapsi) 0;...
        0 0 1];
B = Rb10*([0;e;0]+Rb21*B);

D = [d;0;0];
Rd21 = [cos(deltad2) -sin(deltad2) 0;...
        sin(deltad2) cos(deltad2) 0;...
        0 0 1]; 
Rd20 = [0 1 0;...
        -1 0 0;...
        0 0 1];
D = [0;e;0]+Rd20*([c_d;0;0]+Rd21*D);

F = [f;0;0];
Rf21 = [cos(deltaf2) -sin(deltaf2) 0;...
        sin(deltaf2) cos(deltaf2) 0;...
        0 0 1]; 
Rf20 = [0 1 0;...
        -1 0 0;...
        0 0 1];
F = [0;e;0]+Rf20*([c_f;0;0]+Rf21*F);
%%
daf = F-A;
daf_n =sqrt(daf'*daf);
daf0 = subs(daf_n,[xi1 deltaf2],[0 0]);
daf1 = subs(diff(daf_n,xi1),[xi1 deltaf2],[0 0]);
daf2 = subs(diff(daf_n,deltaf2),[xi1 deltaf2],[0 0]);

k2f = simplifyFraction(-daf1/daf2);
k2fa = subs(k2f,a,0)
k2f_Nb4 = subs(k2f,[deltapsi a c_a c_f f],[2*pi/4 e/2 e/4 e/4 e/2]);
k2f_Nb5 = double(simplifyFraction(subs(k2f,[deltapsi a c_a c_f f],[2*pi/5 e/2 e/4 e/4 e/2])));

dbd = D-B;
dbd_n = sqrt(dbd'*dbd);
dbd0 = subs(dbd_n,[xi3 deltad2],[0 0]);
dbd1 = subs(diff(dbd_n,xi3),[xi3 deltad2],[0 0]);
dbd2 = subs(diff(dbd_n,deltad2),[xi3 deltad2],[0 0]);

k2d = simplifyFraction(-dbd1/dbd2);
solb = solve(k2d==0,b)
solcb = solve(subs(k2d,b,0)==0,c_b)

k2d_Nb4 = subs(k2d,[deltapsi b c_b c_d d],[2*pi/4 e/2 e/4 e/4 e/2]);
k2d_Nb5 = double(simplifyFraction(subs(k2d,[deltapsi b c_b c_d d],[2*pi/5 e/2 e/4 e/4 e/2])));
%%
R = [0 k2d;k2f 0];
C_Ri = R'*[1 -1]'*[1 -1]*R;
syms C_d1 C_d2 C_ed C_d
C_Risym = subs(C_Ri,[C_Ri(1,1) C_Ri(2,2)], [C_d1 C_d2]);
% C_Ri = subs(C_Ri,[a b], [0 0]);
if isequal(C_Ri(2,1),C_Ri(1,2))
    C_Risym = subs(C_Risym,[C_Risym(1,2) C_Risym(2,1)], [C_ed C_ed]);
end
disp(C_Risym)
%%
Nb = 5;
elem = 1:Nb;
C_R = sym(zeros(Nb));
C_Rsym = sym(zeros(Nb));
dpsi = 2*pi/Nb;
for ii=1:Nb
    psi(ii,1) = (ii-1)*dpsi;
    nxt = mod(ii+1,Nb)+1;
    dof(:,ii)=[ii;nxt];
    nxt = mod(ii,Nb)+1;
    C_R(dof(:,ii),dof(:,ii))=C_R(dof(:,ii),dof(:,ii))+C_Ri;
    C_Rsym(dof(:,ii),dof(:,ii))=C_Rsym(dof(:,ii),dof(:,ii))+C_Risym;
    %C_R(dof(:,ant),dof(:,ant))=C_R(dof(:,ant),dof(:,ii))+C_Ri;
end
C_Rsym = subs(C_Rsym,C_d1+C_d2,C_d)
C_R = simplify(C_R,steps=50);

C_R_double = simplifyFraction(subs(C_R,[b c_b a c_a c_d d c_f f],[e/2 e/4 e/2 e/4 e/4 e/2 e/4 e/2]));
C_R_double = double(simplify(subs(C_R_double,[deltapsi],[2*pi/Nb])));

T = rewrite(MBC(Nb),'expandsum');
C_NRsym = simplify(T'*C_Rsym*T);
C_NRsym = 1/Nb*C_NRsym;
disp(C_NRsym)

