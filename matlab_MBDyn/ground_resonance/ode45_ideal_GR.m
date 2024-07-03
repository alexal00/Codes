% close all; clear all; clc
%%
function [sol,t] = ode45_ideal_GR(omega,Nb,y0,mode,file)
% omega = 3;
T = 2*pi/omega;
tspan = [0 10];
% Nb = 4;
airframe = 2;
dof = Nb+airframe;
% IC = 3;
% if IC==0
%     x_0 = 0; x_d0 = 0;
%     y_0 = 0; y_d0 = 0;
%     if Nb ==4
%         xi_R0 = 0*[0.01, -0.004, 0.005, 0.002]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     elseif Nb==5
%         xi_R0 = 0*[0.01, -0.004, 0.005, 0.002, 0]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     end
% elseif IC ==1
%     x_0 = 0.1; x_d0 = x_0/T;
%     y_0 = 0; y_d0 = y_0/T;
%     if Nb ==4
%         xi_R0 = [0.01, -0.004, 0.005, 0.002]';
%         xi_R0d = zeros(Nb,1)+xi_R0/T;
%     elseif Nb==5
%         xi_R0 = [0.01, -0.004, 0.005, 0.002, 0]';
%         xi_R0d = zeros(Nb,1)+xi_R0/T;
%     end
% elseif IC ==2
%     x_0 = 0.1; x_d0 = 0*x_0/T;
%     y_0 = 0; y_d0 = 0*y_0/T;
%     if Nb ==4
%         xi_R0 = [0.01, -0.004, 0.005, 0.002]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     elseif Nb==5
%         xi_R0 = [0.01, -0.004, 0.005, 0.002, 0]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     end
% elseif IC ==3
%     x_0 = 0; x_d0 = 0*x_0/T;
%     y_0 = 0; y_d0 = 0*y_0/T;
%     if Nb ==4
%         xi_R0 = [0.01, -0.004, 0.005, 0.002]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     elseif Nb==5
%         xi_R0 = [0.01, -0.004, 0.005, 0.002, 0]';
%         xi_R0d = zeros(Nb,1)+0*xi_R0/T;
%     end
% end
% y0 = [xi_R0' x_0 y_0 xi_R0d' x_d0 y_d0];
% y0 = [zeros(1,Nb) 0.1 0 zeros(1,Nb) 0 0];
% xi_R0_deg = xi_R0*180/pi;
deltapsi = 2*pi/Nb;
% psi_0 = [deltapsi*(1:Nb)]';
% xi_NR0 = MBC(psi_0,xi_R0);

% Select damping type
% 0 : Blade to hub
% 1 : Interblade
% 2 : Inter-2-blade
% mode = 1;
[t,sol] = ode45(@(t,y) odefun(t,y,omega,mode,Nb),tspan,y0);
% load("ode_results\v4_mode2_Nb4_xi0.mat","t",'sol')
% load("mode2Nb5_sol.mat")
%% Solution in the rotating frame
airframe = 2;
xi = sol(:,1:Nb);
hub = sol(:,Nb+1:Nb+airframe);
% figure(1) 
% for ii = 1:Nb
% plot(t, xi(:,ii),'DisplayName',['$\xi_' num2str(ii) '$']); hold on
% end
% legend('show','interpreter','latex')
% grid on
% 
% %% Solution in the non-rotating frame
% psi_a = zeros(Nb,length(t));
% for jj = 1:Nb
% psi_a(jj,:)=omega*t'+(jj-1)*deltapsi;
% end
% xi_NR_a = MBC(psi_a,xi');
% 
% qNR = NRdof(Nb);
% figure(2)
% subplot(1,2,1)
% % Plot collective DOF
% if mod(Nb,2)==0
%     dof_col = [1 Nb];
%     plot(t, xi_NR_a(1,:),'DisplayName',qNR{1}); hold on
%     plot(t, xi_NR_a(end,:),'DisplayName',qNR{end}); hold on
% else
%     dof_col = 1;
%     plot(t, xi_NR_a(1,:),'DisplayName',qNR{1}); hold on
% end
% xlabel('t [s]','Interpreter','latex')
% ylabel('$q_{NR}$ [rad]','Interpreter','latex')
% legend('show','interpreter','latex')
% grid on
% dof_cyc = setdiff(1:Nb,dof_col);
% subplot(1,2,2)
% for jj = 1:length(dof_cyc)
% plot(t', xi_NR_a(dof_cyc(jj),:),'DisplayName',qNR{dof_cyc(jj)}); hold on
% end
% xlabel('t [s]','Interpreter','latex')
% ylabel('$q_{NR}$ [rad]','Interpreter','latex')
% legend('show','interpreter','latex')
% grid on

save(file,"t","sol")