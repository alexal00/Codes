close all; clear all; clc
%%
set(groot,"defaultLineLineWidth",2,'DefaultAxesTickLabelInterpreter', 'latex',...
    'DefaultTextInterpreter','latex','DefaultLegendInterpreter','latex');
%%
t = linspace(0,3,301);
dt = t(2)-t(1);
Omega = 2*pi;


syms time x_d0 x_0 T b a
v = a-b*heaviside(time-T);
acc = diff(v,time);
x = int(v,time);
% v2 = diff(x,time);
% x_0 = 1;
% x_d0 = 2;
% T = 2*pi/Omega;

v_t = matlabFunction(subs(v,[T a b], [2*pi/Omega (2+1*Omega/(2*pi)) 1*Omega/(2*pi)]));
a_t = matlabFunction(subs(acc,[T a b], [2*pi/Omega (2+1*Omega/(2*pi)) 1*Omega/(2*pi)]));
% v_t2 = matlabFunction(subs(v2,[T a b], [2*pi/Omega (2+1*Omega/(2*pi)) 1*Omega/(2*pi)]));
x_t = matlabFunction(subs(x,[T a b], [2*pi/Omega (2+1*Omega/(2*pi)) 1*Omega/(2*pi)]));
figure
subplot(3,1,1)
plot(t,x_t(t),'DisplayName','Position'); hold on
plot(2*pi/Omega,(2+1*Omega/(2*pi))*2*pi/Omega,'ro','HandleVisibility','off'); hold on
% yline((2+1*Omega/(2*pi))*2*pi/Omega,'r--','HandleVisibility','off');
% xline(2*pi/Omega,'r--','HandleVisibility','off');
text(2*pi/Omega+0.1,(2+1*Omega/(2*pi))*2*pi/Omega-0.8,'($T$,$(\dot{x}_0+\frac{x_0}{T})T$)')
grid on
xlabel('$t$ [s]')
ylabel('$x$ [m]')
% legend('Location','best')

subplot(3,1,2)
plot(t,v_t(t),'DisplayName','Velocity'), hold on
% plot(t,(t<2*pi/Omega),'DisplayName','Step, t$\leq$T');
plot(2*pi/Omega-dt,(2+1*Omega/(2*pi)),'ro','HandleVisibility','off')
plot(2*pi/Omega+dt,2,'ro','HandleVisibility','off')
text(2*pi/Omega+0.1,(2+1*Omega/(2*pi)),'($T-dt$,$\dot{x}_0+\frac{x_0}{T}$)')
text(2*pi/Omega,(2)-0.5,'($T+dt$,$\dot{x}_0$)')
xlabel('$t$ [s]')
ylabel('$v_x$ [m/s]')
ylim([0 (2+1*Omega/(2*pi))+0.5])
% legend('Location','best')
grid on

subplot(3,1,3)
a_vec = a_t(t);
idx = a_vec == -Inf; % find Inf
a_vec(idx) = -1*Omega/(2*pi);     % set Inf to finite value
stem(t,a_vec,'HandleVisibility','off'); hold on
plot(2*pi/Omega,-1*Omega/(2*pi),'ro','HandleVisibility','off'); hold on
text((2*pi)/Omega+0.1,-1*Omega/(2*pi)+0.2,'($T$,$-x_0/T$)')
xlabel('$t$ [s]')
ylabel('$a_x$ [m/s\textsuperscript{2}]')
ylim([-1*Omega/(2*pi)-0.2 0])
grid on