% Author: Alejandro Alvaro, 2023-2024
% Obtain geometrical parameters to include in the generalized Deutsch criterion to size
% the lead-lag dampers for the different options investigated
close all; clear all; clc
%%
e = 0.3048;     % [m], blade hinge offset from Hammond [1974]
% Geometrical consideration
Nb = [4,5,6];
opt = [1 2 3];
for nb = Nb
    disp(['Number of blades ' int2str(nb)])
    for ii = 1 : length(opt)
        if opt(ii) ==1 % original case
            a = e/2;
            b = e/2;
            ca = e/4;
            cb = e/4;
            d = e/2;
            f = e/2;
            cd = e/4;
            cf = e/4;
        elseif opt(ii)==2 % modified case
            a = e/4;
            b = e/4;
            ca = e/4;
            cb = e/4;
            d = e/4;
            f = e/4;
            cd = e/2;
            cf = e/2;
         elseif opt(ii)==3 % modified case
            a = e/4;
            b = e/4;
            ca = e/4;
            cb = e/4;
            d = e/4;
            f = e/4;
            cd = e/4;
            cf = e/4;
        end
        disp(['Option number, opt=' int2str(opt(ii))])
        disp('The geometrical constants are:')
        disp('i2b:')
        [Kf,Kd] = i2bfun(e,a,b,ca,cb,d,f,cd,cf,nb);
        disp(['Kxif =' num2str(Kf)])
        disp(['Kxid =' num2str(Kd)])
        disp('ib:')
        [Ka,Kb] = ibfun(e,a,b,ca,cb,nb);
        disp(['KxiA =' num2str(Ka)])
        disp(['KxiB =' num2str(Kb)])
    end
end