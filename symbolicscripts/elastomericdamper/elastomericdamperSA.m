close all; clear all; clc
%%
setPlot
%%
figure("Name",'tanh')
plot([-10:1/1000:10],sign([-10:1/1000:10]),'DisplayName','$\mathrm{sgn}(x)$')
alpha_v = logspace(-2,2,5); hold on
for ii = 1 : length(alpha_v)
    plot([-10:1/1000:10],tanh(alpha_v(ii).*[-10:1/1000:10]),'DisplayName',['$\tanh{(' num2str(alpha_v(ii)) 'x)}$'])
end
xlabel('$x$')
ylabel('$y$')
grid on
legend

%%
Ktheta = 14000;
fd = 4067.5;
alpha = 10;
Omega = 40;
T = 2*pi/Omega;
dt = T/256;
t = 0:T/256:T;
A = 1*pi/180;
theta = A*sin(Omega*t);
thetap = Omega*A*cos(Omega*t);
Fd = Ktheta*theta + fd*tanh(alpha*thetap);
Fdsgn = Ktheta*theta + fd*sign(alpha*thetap);
figure('Name','CLtanh')
subplot(1,2,1)
plot(theta*180/pi,Fd); hold on
% plot(theta*180/pi,Fdsgn);

xlabel('$\phi$ [deg]')
ylabel('$F_d$ [Nm]')
grid on
subplot(1,2,2)
plot(thetap*180/pi,Fd);hold on
% plot(thetap*180/pi,Fdsgn);
xlabel('$\dot{\phi}$ [deg/s]')
ylabel('$F_d$ [Nm]')
grid on
%% Sensitivy analysis of Ktheta on CL
Kthetavec = [0 logspace(1,5,5)];
figure('Name','CLvsKtheta')
tiledlayout(1,2)
ax1 = nexttile;
xlabel(ax1,'$\phi$ [deg]')
ylabel(ax1,'$F_d$ [Nm]')
grid(ax1,"on")
hold(ax1,"on")
% legend(ax1)
ax2 = nexttile;
xlabel(ax2,'$\dot{\phi}$ [deg/s]')
ylabel(ax2,'$F_d$ [Nm]')
grid(ax2,"on")
hold(ax2,"on")
% legend(ax2,'Location','bestoutside')
lg  = legend(ax2,'Orientation','Vertical','NumColumns',5); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout
for ii = 1 :length(Kthetavec)
    Fd = Kthetavec(ii)*theta + fd*tanh(alpha*thetap);
    plot(ax1,theta,Fd,'DisplayName',['$K_{\phi}$= ' num2str(Kthetavec(ii))]);
    plot(ax2,thetap,Fd,'DisplayName',['$K_{\phi}$= ' num2str(Kthetavec(ii))]);
end
%% Sensitivity analysis of f and alpha on Ceq

omega = 2*pi;
t = linspace(0,2*pi/omega,201);
A = 10*pi/180;
alpha_v = linspace(1,300,101);
figure('Name','Ceqvsalpha')
f_cte = 1000;
C_eq = zeros(1,length(alpha_v));
for ii=1:length(alpha_v)
    C_eq(ii) = calculate_C_eq_given_f(f_cte,omega,A,alpha_v(ii),t);
end
C_eqsgn = calculate_C_eq_given_f_sgn(f_cte,omega,A,t);

plot(alpha_v,C_eq,'DisplayName','$\bar{f}\tanh{(\alpha\dot{x})}$'); hold on
plot(alpha_v,C_eqsgn*ones(size(alpha_v)),'DisplayName','$\bar{f}\mathrm{sgn}(\dot{x})$');
xlabel('$\alpha$')
ylabel('$C_{eq}$')
legend
grid on

figure('Name','Ceqvsf')
alpha_cte = 50;
f_v = linspace(250,5000,51);
C_eq = zeros(1,length(f_v));
for ii=1:length(f_v)
    C_eq(ii) = calculate_C_eq_given_f(f_v(ii),omega,A,alpha_cte,t);
end
dCeqdf = (C_eq(2)-C_eq(1))/(f_v(2)-f_v(1));
Kappa = dCeqdf/(A*omega);
plot(f_v,C_eq,'DisplayName','$\bar{f}\tanh{(\alpha\dot{x})}$'); hold on
% plot(500,calculate_C_eq_given_f(500,omega,A,alpha_cte,t),'Marker','o')
xlabel('$\bar{f}$')
ylabel('$C_{eq}$')
% legend
grid on

return
%% Deutsch criterion for friction damper
syms f alpha x
Fd = f*tanh(alpha*x);
Fd0 = subs(Fd,x,0);
Fd1 = subs(diff(Fd,x),x,0);

%% Auxiliary functions
% Define the function to calculate C_eq given a value of f
function C_eq = calculate_C_eq_given_f(f, omega, A, alpha, t)
    thetap = omega * A * cos(omega * t);
    F = @(nu) f * tanh(alpha * nu);
    Ed = trapz(t, F(thetap) .* thetap);
    C_eq = Ed / (pi * omega * A^2);
end

function C_eq = calculate_C_eq_given_f_sgn(f, omega, A, t)
    thetap = omega * A * cos(omega * t);
    F = @(nu) f * sign(nu);
    Ed = trapz(t, F(thetap) .* thetap);
    C_eq = Ed / (pi * omega * A^2);
end