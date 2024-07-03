% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model for GR and visualize output
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;            % Figure id

%% Helicopter data and envrionment variables
%======================= GENERAL HELICOPTER parameters=================
Nb = 5;             % [-], number of blades
omega = 3.;         % [rad/s], rotor angular speed

% Execute or no MBDyn
% 1 : execute MBDyn
% else: do not execute
exec = 0.;

% Geometrical consideration
opt = 1;

e = 0.3048;

ic = 3;
%% MBDyn pathfiles
% path files for MBDyn
pref_in = '/home/ale/thesis/GR_MBDyn/';
pref_out = '/mnt/c/Users/aasal/Documents/Universidad/Master/Thesis/Codes/matlab_MBDyn/ground_resonance/';
sf_in = {'alvaro_model/v10/'};
% sf_in = {'cassoni_model/'};
% filename = 'rotor_aero_i2b_Nb4.mbd';
damping = {'std' 'ib' 'i2b'};
% damping = {'i2b'};
% damping = {'ib'};
sf_out = {'MBDynv10/'};
%% MBDyn variables
geom = calculateKgeom(opt,e,damping,Nb);

variables =['"' 'const real OMEGA_100 = ' num2str(omega) ...
    '; const real A_BAR = ' num2str(geom.a) ...
    '; const real B_BAR = ' num2str(geom.b) ...
    '; const real C_BAR = ' num2str(geom.cd) ...
    '; const real C_BARA = ' num2str(geom.ca) ...
    '; const real C_BARB = ' num2str(geom.cb) ...
    '; const real D_BAR = ' num2str(geom.d) ...
    '; const real F_BAR = ' num2str(geom.f) ...
    '; const real GEOM = ' num2str(geom.Kxidelta^2) ...
    '; const real GEOMIB = ' num2str(geom.Kxil^2) ...
    '"'];
% setenv('MBDYNVARS', variables);
%% Execute MBDyn

for i = 1:length(damping)
    filename = ['rotor_' damping{i} '_Nb' int2str(Nb) '_IC'  num2str(ic) '.mbd'];
    % filename = 'GR_rotor.mbd';
    outfile = [damping{i} '_Nb' int2str(Nb) '_Hammond_GR'];
    % Enviroment variables for MBDyn
    % Select damping type
    % 0 : Blade to hub
    % 1 : Interblade
    % 2 : Inter-2-blade
    if contains(filename,'i2b')
        mode = 2.;
    elseif contains(filename,'ib')
        mode = 1;
    else
        mode = 0;
    end

    % ODE results for comparisson
    file = ['./GRode_results/v11_mode' int2str(mode) '_Nb' int2str(Nb) '_IC' num2str(ic) '_O3.mat'];


    for kk = 1 : length(sf_in)
    % Input file path
    in = [pref_in sf_in{kk}];
    % Output file path seen from wsl pov
    out = [pref_out sf_out{kk}];
    % Output file path seen from the windows pov
    outwd = sf_out{kk};
    % Execute MBDYN
    [params, mbdynsol] = executembdyn(in,out,outwd,filename,outfile,fid,variables,exec);
    outpath = [sf_out{kk} outfile '.nc'];
    % return
    % continue

    % Nb = params.nb;         % Number of blades
    airframe = 2;           % Airframe dof

    % MBC transformation

    % Preallocate matrices
    xi_R_MB = zeros(4,length(mbdynsol.t));
    psi = zeros(4,length(mbdynsol.t));
    for jj = 1: Nb
        xi_R_MB(jj,:) =  mbdynsol.(['xi_' num2str(jj)])';
        psi(jj,:) = mbdynsol.(['psi_' num2str(jj)])';
    end
    xi_NR_MB = MBC(psi,xi_R_MB);

    % Transfrom the output structure to table
    mbdyn_tab = struct2table(mbdynsol);

    % Check relationship between lead-lag angle and damper deformation
    % Only valid for the i2b damper
    row = floor(Nb/2);
    col = ceil(Nb/2);
    aux = params.index;
    % if mode==2
        % fhandle = figure(fid); fid = fid+1;
        % fhandle.Name = ['Nb' int2str(Nb) damping{i} 'xidelta'];
        % for ii=1:Nb
        %     ant = mod(ii - 2, Nb) + 1;  % Calculate previous index
        %     nxt = mod(ii, Nb) + 1;  % Calculate next index
        %     nm1=['xi_' num2str(ii)];
        %     leg1 = latexfmt(nm1);
        %     xi = mbdyn_tab.(nm1);
        %     % Arm 1 of damper i-1
        %     nm2=['delta_D' num2str(ant)];
        %     leg2 = latexfmt(nm2);
        %     delta_D = mbdyn_tab.(nm2);
        % 
        %     % Arm 2 of damper i+1
        %     nm3=['delta_F' num2str(nxt)];
        %     leg3 = latexfmt(nm3);
        %     delta_F = mbdyn_tab.(nm3);
        % 
        %     % LS fit for F and D
        %     % k_F(ii) = constant(xi,delta_F);
        %     % k_D(ii) = constant(xi,delta_D);
        %     % text = {['$\xi_i=k_F\delta_{Fi+1}$, $k_F$=' num2str(k_F(ii))],['$\xi_i=k_{D}\delta_{Di-1}$, $k_D$=' num2str(k_D(ii))]};
        %     % ax = subplot(row,col,ii)
        %     % functionplot(ax,mbdyn_tab.t,xi,step,Color{ii},Marker{ii},'none',leg1)
        %     % functionplot(ax,mbdyn_tab.t,delta_F,step,Color{nxt},Marker{nxt},'none',leg2)
        %     % functionplot(ax,mbdyn_tab.t,delta_D,step,Color{ant},Marker{ant},'none',leg3)
        %     % % gtext(text,'interpreter','latex','HorizontalAlignment', 'center', 'FontSize', 12);
        %     % % plot(mbdyn_tab.t,mbdyn_tab.(nm1),'DisplayName',leg1); hold on
        %     % % plot(mbdyn_tab.t,mbdyn_tab.(nm2),'DisplayName',leg2); hold on
        %     % xlabel('t [s]','Interpreter','latex')
        %     % ylabel('$\xi$, $\delta$ [rad]','Interpreter','latex')
        %     % grid on
        %     % legend('show','interpreter','latex')
        % end
        % fhandle = figure(fid); fid = fid+1;
        % fhandle.Name = ['Nb' int2str(Nb) damping{i} 'Deltadelta'];
        % for ii=1:Nb
        %     ant = mod(ii - 2, Nb) + 1;  % Calculate previous index
        %     nxt = mod(ii, Nb) + 1;  % Calculate next index
        %     % Relative orientation between both
        %     nm1=['Deltaxi_' num2str(ii)];
        %     leg1 = ['$\Delta\delta_' num2str(ii) '$'];
        %     dxi = mbdyn_tab.(nm1);
        %     % Arm 1 of the damper
        %     nm2=['delta_D' num2str(ii)];
        %     leg2 = latexfmt(nm2);
        %     delta_D = mbdyn_tab.(nm2);
        %     % Arm 2 of the damper
        %     nm3=['delta_F' num2str(ii)];
        %     leg3 = latexfmt(nm3);
        %     delta_F = mbdyn_tab.(nm3);
        % 
        %     ax = subplot(row,col,ii)
        %     % functionplot(ax,mbdyn_tab.t,dxi,step,Color{ii},Marker{ii},'none',leg1)
        %     functionplot(ax,mbdyn_tab.t,delta_D,step,Color{ant},Marker{ant},'none',leg2)
        %     functionplot(ax,mbdyn_tab.t,delta_F,step,Color{nxt},Marker{nxt},'none',leg3)
        %     xlabel('t [s]','Interpreter','latex')
        %     ylabel('$\delta$ [rad]','Interpreter','latex')
        %     grid on
        %     legend('show','interpreter','latex')
        % end
    % elseif mode==0
    %     fhandle = figure(fid); fid = fid+1;
    %     fhandle.Name = ['Nb' int2str(Nb) damping{i} 'damper'];
    %     for ii=1:Nb
    %         colororder({'k','k'})
    %         label = int2str(10000 + 1000*ii+30+1);
    %         xi = ncread(outpath,['elem.joint.' label '.E'],[3 1],[1 inf]);
    %         xi = xi(aux:end);
    %         xi_d = ncread(outpath,['elem.joint.' label '.Omega'],[3 1],[1 inf]);
    %         xi_d = xi_d(aux:end);
    %         ax = subplot(row,col,ii)
    %         yyaxis left
    %         functionplot(ax,mbdyn_tab.t,xi,step,Color{ii},Marker{ii},'none',['$\xi_' num2str(ii) '$'])
    %         ylabel('$\xi$ [rad]')
    %         yyaxis right
    %         functionplot(ax,mbdyn_tab.t,xi_d,step,Color{ii+1},Marker{ii},'none',['$\dot{\xi}_' num2str(ii) '$'])
    %         % plot(mbdyn_tab.t,xi_d,'DisplayName',['$\dot{\xi}_' num2str(ii) '$']); hold on
    %         ylabel('$\dot{\xi}$ [rad/s]')
    %         xlabel('t [s]')
    %         legend
    %         grid on
    % 
    %     end
    % end
    %% ODE solution comparisson
    % aux = find(mbdyn_tab.t>=0,1);
    v_x = ncread(outpath,['node.struct.' num2str(100) '.XP']);
    vx_0 = v_x(1,aux);
    v_y = ncread(outpath,['node.struct.' num2str(200) '.XP']);
    vy_0 = v_y(2,aux);
    % ODE parameters
    tspan = [0 2];
    
    for ii=1:Nb
        actual = int2str(10000+1000*ii);
        xi_0(ii) = mbdyn_tab.(['xi_' num2str(ii)])(1);
        xip(ii,:) = ncread(outpath, ['node.struct.',actual,'.Omega'],[3, 1], [1, inf])-omega;
        xip_0(ii) = xip(ii,aux-1);
    end
    x_0 = mbdyn_tab.hub_x(1);
    y_0 = mbdyn_tab.hub_y(1);
    % Initial conditions
    y0 = [xi_0 x_0 y_0 xip_0 vx_0 vy_0];
    disp('The initial conditions are:')
    disp(y0)
    if exist(file,"file")
        load(file,"t","sol")
    else
        [sol,t] = ode45_ideal_GR(omega,Nb,y0,mode,file);
    end
    % fhandle=figure(fid); fid=fid+1;
    % fhandle.Name = ['Nb' int2str(Nb) damping{i} 'IC'];
    % subplot(2,1,1)
    % plot(mbdyn_tab.t(1:aux-1),v_x(1,1:aux-1),'DisplayName','$v_x$ [m/s]'); hold on
    % plot(mbdyn_tab.t(1:aux-1),v_y(2,1:aux-1),'DisplayName','$v_y$ [m/s]'); hold on
    % xlabel('t [s]','Interpreter','latex')
    % ylabel('$v_x$, $v_y$ [m/s]','Interpreter','latex')
    % legend('show','interpreter','latex')
    % grid on
    % subplot(2,1,2)
    % for ii=1:Nb
    %     nm = ['$\dot{\xi}_' num2str(ii) '$'];
    %     plot(mbdyn_tab.t(1:aux),xip(ii,1:aux),'DisplayName',nm); hold on
    % end
    % xlabel('t [s]','Interpreter','latex')
    % ylabel('$\dot{\xi}$ [rad/s]','Interpreter','latex')
    % legend('show','interpreter','latex')
    % grid on
    
    % Data
    xi = sol(:,1:Nb);
    hub = sol(:,Nb+1:Nb+airframe);
    deltapsi = 2*pi/Nb;
    psi_a = zeros(Nb,length(t));
    for jj = 1:Nb
    psi_a(jj,:)=omega*t'+(jj)*deltapsi;
    end
    xi_NR_a = MBC(psi_a,xi');
    
    
    %% Figures
    % fhandle=figure(fid); fid = fid+1;
    % fhandle.Name = ['Nb' int2str(Nb) damping{i} 'hub'];
    % ax = gca;
    % functionplot(ax,mbdyn_tab.t,mbdyn_tab.hub_x,step,Color{1},Marker{1},'none','hub-x, MBDyn')
    % functionplot(ax,t,hub(:,1),step,Color{1},'none',LineStyle{1},'hub-x, ODE')
    % functionplot(ax,mbdyn_tab.t,mbdyn_tab.hub_y,step,Color{2},Marker{2},'none','hub-y, MBDyn')
    % functionplot(ax,t,hub(:,2),step,Color{2},'none',LineStyle{2},'hub-y, ODE')
    % xlabel('t [s]','Interpreter','latex')
    % ylabel('x, y [m]','Interpreter','latex')
    % legend('show','interpreter','latex')
    % grid on

    fhandle=figure(fid); fid = fid+1; 
    fhandle.Name = ['Nb' int2str(Nb) damping{i} 'qRIC' num2str(ic)];
    ax = gca;
    for jj = 1:Nb
        functionplot(ax,mbdyn_tab.t,mbdyn_tab.(['xi_' num2str(jj)]),step,Color{jj},Marker{jj},'none',['$\xi_' num2str(jj) '$, MBDyn'])
        plot(ax,t, xi(:,jj),'DisplayName',['$\xi_' num2str(jj) '$'],'LineStyle',LineStyle{jj},'Color',Color{jj}); hold on
    end
    xlabel('t [s]','Interpreter','latex')
    ylabel('$\xi$ [rad]','Interpreter','latex')
    legend('Location','best')
    grid on
    
    fhandle=figure(fid); fid = fid+1;
    fhandle.Name = ['Nb' int2str(Nb) damping{i} 'qNRIC' num2str(ic)];
    qNR = NRdof(Nb);
    ax1 = subplot(1,2,1);
    % Plot collective DOF
    qNRmax = max(xi_NR_MB(:,:),[],"all");
    qNRmin = min(xi_NR_MB(:,:),[],"all");
    if mod(Nb,2)==0
        dof_col = [1 Nb];
        functionplot(ax1,mbdyn_tab.t,xi_NR_MB(1,:),step,Color{1},Marker{1},'none',[qNR{1} ', MBDyn'])
        plot(ax1,t, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
        functionplot(ax1,mbdyn_tab.t,xi_NR_MB(end,:),step,Color{2},Marker{2},'none',[qNR{end} ', MBDyn'])
        plot(ax1,t, xi_NR_a(end,:),'DisplayName',qNR{end},'Color',Color{2}); hold on
    else
        dof_col = 1;
        functionplot(ax1,mbdyn_tab.t,xi_NR_MB(1,:),step,Color{1},Marker{1},'none',[qNR{1} ', MBDyn'])
        plot(ax1,t, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
    end
    xlabel('t [s]','Interpreter','latex')
    ylabel('$q_{NR}$ [rad]','Interpreter','latex')
    legend('Location','best')
    ylim([qNRmin+0.5*qNRmin qNRmax+0.5*qNRmax])
    grid on

    dof_cyc = setdiff(1:Nb,dof_col);
    ax2 = subplot(1,2,2);
    for jj = 1:length(dof_cyc)
        functionplot(ax2,mbdyn_tab.t,xi_NR_MB(dof_cyc(jj),:),step,Color{jj},Marker{jj},'none',[qNR{dof_cyc(jj)} ', MBDyn'])
        plot(ax2,t', xi_NR_a(dof_cyc(jj),:),'DisplayName',qNR{dof_cyc(jj)},'Color',Color{jj},'LineStyle',LineStyle{jj}); hold on
    end
    xlabel('t [s]','Interpreter','latex')
    ylabel('$q_{NR}$ [rad]','Interpreter','latex')
    legend('Location','best')
    grid on
    ylim([qNRmin+0.5*qNRmin qNRmax+0.5*qNRmax])
    
    fid = fid+10;
    % clear t sol
    end
end
%%
savefigures