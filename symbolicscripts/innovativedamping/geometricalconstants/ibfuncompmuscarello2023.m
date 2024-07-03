% Author: Alejandro Alvaro, 2023-2024
% Comparisson between equivalent formulations to obtain the geometrical
% constant for the inter-blade arrangement.
% 1: Proposed method from Thesis
% 2: Extended version from Muscarello & Quaranta [2023]
close all; clc; clear all
%%
xiE = 0;
R = 5.5;
Nb = 4;
e = 0.04*R;

h = [0.5 1 2];
o = linspace(0,1,21);
dpsi = 2*pi/Nb;
phi = (pi-dpsi)/2;
Le=  2*e*sin(dpsi/2);
for ii=1:length(h)
    for jj=1:length(o)
        cb = e.*h(ii).*(1+o(jj));
        ca = e.*h(ii).*(1-o(jj));
        
        gamma = atan(-(ca*sin(xiE-phi)+cb*sin(xiE+phi))/((Le+ca*cos(xiE-phi)+cb*cos(xiE+phi))));
        Cxid = ca^2*sin(xiE-phi+gamma)^2+cb^2*sin(xiE+phi+gamma)^2;
        Cxied = ca*cb*sin(xiE-phi+gamma)*sin(xiE+phi+gamma);
        [~,~,C_R] = ibfun(e,0,0,ca,cb,Nb);
        Cxidibfun = C_R(1,1);
        Cxiedibfun = C_R(1,2);
        if abs(Cxid-Cxidibfun)<=1e-5 && abs(Cxied-Cxiedibfun)<=1e-5
            disp('Correct formulation')
        else
            disp('Incorrect formulation')
        end
    end
end