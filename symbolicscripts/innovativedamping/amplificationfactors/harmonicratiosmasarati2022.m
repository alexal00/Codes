% Author: Alejandro Alvaro, 2023-2024
% Obtain analytical expression for the harmonic amplitude expected in the
% ideal case where only lead-lag motions are allowed and no kinematic
% ocuplings with additional dof are present.
close all; clear all; clc
%%
syms n dpsi real
syms Nb integer

ib_harm = rewrite((exp(1i*n*dpsi)-1),"sincos");
ib_harm = expand(simplify(sqrt(ib_harm*conj(ib_harm))))
ib_harm = subs(ib_harm,dpsi,2*pi/Nb);
i2b_harm = rewrite((exp(1i*2*n*dpsi)-1),"sincos");
i2b_harm = simplify(sqrt(i2b_harm*conj(i2b_harm)))
i2b_harm = subs(i2b_harm,dpsi,2*pi/Nb);

blades = [4,5,6];
for ii=1:length(blades)
    ib_harmfun = matlabFunction(ib_harm);
    i2b_harmfun = matlabFunction(i2b_harm);
    
    n = 0:7;

    figure
    bar(n,ib_harmfun(blades(ii),n),'DisplayName','ib harmonics'); hold on
    bar(n,i2b_harmfun(blades(ii),n),'DisplayName','i2b harmonics'); hold on
    xlabel('Harmonic')
    ylabel('Amplif')
    legend
end
%% Individual plots
% Obtain individual bar plots for each of the arrangements
Nb = 5;
dpsi = 2*pi/Nb;
n = 1:7;

fact1 = sqrt(2*(1-cos(2*n*dpsi)));
fact2 = fact1/(2*(1-cos(2*dpsi)));

figure('Name','i2bharmonics')
subplot(1,2,1)

bar(n,fact1); hold on
title('$\Xi^{i2b}=\sqrt{2(1-\cos{(2n\Delta\psi)})}$','Interpreter','latex')
subplot(1,2,2)
bar(n,fact2); hold on
title('$\mathcal{M}^{i2b}=\frac{\sqrt{2(1-\cos{(2n\Delta\psi)})}}{2(1-\cos{(2\Delta\psi)})}$','Interpreter','latex')


fact1 = sqrt(2*(1-cos(n*dpsi)));
fact2 = fact1/(2*(1-cos(dpsi)));
figure('Name','ibharmonics')
subplot(1,2,1)

bar(n,fact1); hold on
title('$\Xi^{ib}=\sqrt{2(1-\cos{(n\Delta\psi)})}$','Interpreter','latex')
subplot(1,2,2)
bar(n,fact2); hold on
title('$\mathcal{M}^{ib}=\frac{\sqrt{2(1-\cos{(n\Delta\psi)})}}{2(1-\cos{(\Delta\psi)})}$','Interpreter','latex')