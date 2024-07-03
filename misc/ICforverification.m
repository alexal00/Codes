% Author: Alejandro Alvaro, 2023-2024
% obtain initial conditions detailed in Table D.1 from the thesis
% ¨A Comparative Study of the Impact of Innovative Lead-lag Damping 
% Configurations on Helicopter Rotor Stability and Loads¨
close all; clear all; clc
%% Initial conditions
Nb = 4;
dpsi = 2*pi/4;

psi_0=[dpsi*(1:Nb)]';

% Only collective
xi_NR0=[0.01;0;0;0];
xi_R0 = MBC_inv(psi_0,xi_NR0);
disp('Collective only IC')
disp(xi_R0);

% Only reactionless
xi_NR0=[0;0;0;0.01];
xi_R0 = MBC_inv(psi_0,xi_NR0);
disp('Reactionless only IC')
disp(xi_R0);

% Only cosine cyclic
xi_NR0=[0;0.01;0;0];
xi_R0 = MBC_inv(psi_0,xi_NR0);
disp('Cyclic, xi_1c, only IC')
disp(xi_R0);

% Only sine cyclic
xi_NR0=[0;0;0.01;0];
xi_R0 = MBC_inv(psi_0,xi_NR0);
disp('Cyclic, xi_1s, only IC')
disp(xi_R0);

% Only cyclic
xi_NR0=[0;0.02;0.01;0];
xi_R0 = MBC_inv(psi_0,xi_NR0);
disp('Cyclic, xi_1s, only IC')
disp(xi_R0);