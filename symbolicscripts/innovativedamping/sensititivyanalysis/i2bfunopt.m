close all; clear all; clc
%%
setPlot
fid = 1;
%%
R = 5.5;
% e = 0.04*R;
% blades = [4 5 6];
e = 0.04*R;
Nb = [4,5,6];
opt = [1 2 3];
labels = ["Option 1" "Option 2" "Option 3"];
fh = figure(fid); fid=fid+1;
fh.Name = 'DampEffi2b';
ax = gca;
hold(ax,"on");
idx = 1;
for nb = Nb
    
    dpsi = 2*pi/nb;
    deutsch = 2*(1-cos(2*dpsi));
    
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
        [Kf,Kd] = i2bfun(e,a,b,ca,cb,d,f,cd,cf,nb);
        E_c(ii,idx) = (Kd^2+Kf^2-2*Kd*Kf*cos(2*dpsi))/deutsch;
    end    
    idx = idx+1;
    
end
bar(ax,labels,E_c);
ylabel('$E_{fc}=C_{\xi}^{i2b}/C_{\xi}^{i2b}|_{ref}$')
legend(ax,{'$N_b$= 4' '$N_b$=5' '$N_b$=6'},'Location','northwest')
grid(ax,"on")
%%
savefigures