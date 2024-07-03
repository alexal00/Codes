% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
% NOTE: the corresponding model is thougth as a ideal situation WITHOUT
% aerodynamics to see if the evaluation of the possible pitch-lag couplings
% is correct.
% This is done making use of the function kinconsti2b.m
%
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;
lastfid = fid;
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%================= FLIGHT CONDITIONS ============================
% % set MBDyn input for hover
V_INF = 0.; % m/s
ALTITUDE = 0.; % m
% NOTE: for convenience, mount the rotor with vertical shaft and put gamma + tau in alpha_H
ALPHA_H = 0.*pi/180; % radian
INPUT_THETA_0 = 10.*pi/180; % radian
INPUT_A_1_H = 0*pi/180; % radian
INPUT_B_1_H = 5.*pi/180; % radian

%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades

omega = 10.;     %[rad/s], rotor angular speed
nper = 60;      % [-], number of periods to run
R = 5.5;        % [m], rotor radius
e = 0.3048;     % [m], hinge position from hub
nuxi = sqrt(3/2*e/(R-e));

% Important hinge numbers for later post-processing
% ------------ NOVEL BLADE MODEL --------------- 
%	X_BLADE_FLAP_OFFSET 
%  <------------> 
% 
%	c_bar	c 
%  <-----|------> 
%  	   d 
%	  <--> 
%				o							^ 
%	   o		|	___________________		|	a 
%		\		|	|					\	| 
%  o-----o======o---|					|	- 
%		/		|	|					|	| 
%	   o		|	--------------------	|	b 
%				0							v 
%	  <--> 
%	   f 
% 
D_NOD = 80;                 % Arm D of angular damper
F_NOD = 90;                 % Arm F of angular damper
HUB = 10000;                % Hub node
BLADE_FLAP_OFFSET = 30;     % Blade hinge

%% Parametric definition of the model
% Activate or deactivate additional blade-hub damper in blades
% 0 : inactive
% 1 : active
% C_b2h = b2hdamp*C_xi(gamma~0.8)
b2hdamp = 1.;

% Swashplate active or not:
% 1: Active swashplate
% else: No-swashplate and pitch constrained
sw = 1;

% Activate or deactivate blade degrees of freedom
% Since the degrees of freedom are given by means of a total joint in MBDyn
% 0 means that the degrre of freedom is not constrained and viceversa for
% 1.
% 0: allowed
% 1: clamped
if sw ==1
    pitch = 0;
else
    pitch = 1;
end
flap = 0;
lag = 0;


% Additional figure naming convention depending on the studied DOF
dof = '';
if pitch ==0
    dof = [dof 'p'];
end
if flap ==0
    dof = [dof 'f'];
end
if lag ==0
    dof = [dof 'l'];
end

% Execute or no MBDyn
% 1 : execute MBDyn
% else: do not execute
exec = 0.;

% Geometrical consideration
opt = 1;

% CL consideration
% 'lve': Linear visco-elastic damper
% 'fr': ¨friction" damper
cl = 'lve';

% Select the casefiles to be studied
damping = {'i2b'};


% damping = {'std' 'ib' 'i2b'};
% damping = {'std'};
% damping = {'ib'};
% damping = {'std' 'ib' 'i2b'};
% damping = {'ib'};
% damping = {'std' 'ib'};

%% MBDyn variables
% ------------Positioning of damper arms---------------------%
geom = calculateKgeom(1,e,damping,Nb);

%% Global definition of variables
variables =['"' 'const real R = ' num2str(R) ...
    '; const real OMEGA_100 = ' num2str(omega) ...
    '; const real N_PER = ' num2str(nper) ...
    '; const real V_INF = ' num2str(V_INF) ...
    '; const real ALTITUDE = ' num2str(ALTITUDE) ...
    '; const real ALPHA_H = ' num2str(ALPHA_H) ...
    '; const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0) ...
    '; const real INPUT_A_1_H = ' num2str(INPUT_A_1_H) ...
    '; const real INPUT_B_1_H = ' num2str(INPUT_B_1_H) ...
    '; const real ACTIVE = ' num2str(b2hdamp) ...
    '; const integer PITCH = ' num2str(pitch) ...
    '; const integer FLAP = ' num2str(flap) ...
    '; const integer LAG = ' num2str(lag) ...
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
setenv('MBDYNVARS', variables);
%% MBDyn pathfiles
% path files for MBDyn
pref_in = '/home/ale/thesis/KINCUP/';
pref_out = '/mnt/c/Users/aasal/Documents/Universidad/Master/Thesis/Codes/kinematic_coupling/';
sf_in = {'v1/'};
model = 'rotor_';
sf_out = {'v1/'};
%% Execute MBDyn
for i = 1:length(damping)
    fignames = [damping{i} int2str(Nb)];
    if sw ==1
        swash = 'sw';
            
    else
        swash = 'nsw';
    end
    filename =  [model damping{i} '_Nb' int2str(Nb) '_' swash  cl '.mbd'];
    disp(['Input file: ' filename])
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


    % Transfrom the output structure to table
    mbdyn_tab = struct2table(mbdynsol);

    % Check relationship between lead-lag angle and damper deformation
    % Only valid for the i2b damper
    if mode==2
        [Ktt,Kbb,Ktb,Kxx] = kinconsti2b(e,geom.a,geom.ca,geom.f,geom.cf,Nb);
        fhandle = figure(fid); 
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'xidelta'];
        tlo(fid) = tiledlayout('flow');fid = fid+1;

        fhandle = figure(fid); 
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'thetadelta'];
        tlo(fid) = tiledlayout('flow'); fid = fid+1;

        fhandle = figure(fid); 
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'thetabeta'];
        tlo(fid) = tiledlayout('flow'); fid = fid+1;

        fhandle = figure(fid);
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'errorKINCOUP'];
        tlo(fid) = tiledlayout('flow');
        for ii=1:Nb
            % Loop indexes for blade numberings
            ant = mod(ii - 2, Nb) + 1;  % Calculate previous index
            nxt = mod(ii, Nb) + 1;  % Calculate next index
            % node_label = int2str(BLADE(ii));
            % Obtain relevant hinge number
            prev = int2str(HUB+1000*ant + D_NOD);               % D arm
            actual = int2str(HUB+1000*ii + BLADE_FLAP_OFFSET);  % blade hinge
            next = int2str(HUB+1000*nxt + F_NOD);               % F arm

            

            % Legend labels
            nm1=['xi_' num2str(ii)];
            leg1 = latexfmt(nm1);
            xi_i = mbdyn_tab.(nm1);
            % Arm 1 of damper i-1
            nm2=['delta_D' num2str(ant)];
            leg2 = latexfmt(nm2);
            delta_D = mbdyn_tab.(nm2);

            nmt =['theta_' num2str(ii)];
            legt = latexfmt(nmt);
            % Arm 2 of damper i+1
            nm3=['delta_F' num2str(nxt)];
            leg3 = latexfmt(nm3);
            delta_F = mbdyn_tab.(nm3);
            
            nmb =['beta_' num2str(ii)];
            legb = latexfmt(nmb);

            % Orientation of damper arms
            % * connection with previous blade
            delta_d = ncread(outpath, ['elem.joint.',prev,'.Phi'],[3, 1], [1, Inf]);
            % * connection with next blade
            delta_f = ncread(outpath, ['elem.joint.',next,'.Phi'],[3, 1], [1, Inf]);
            % Lead lag angle
            xii = ncread(outpath, ['elem.joint.',actual,'.Phi'],[3, 1], [1, Inf]);
            xiip = ncread(outpath, ['elem.joint.',actual,'.Omega'],[3, 1], [1, Inf]);
            % Pitch angle and velocity
            theta = ncread(outpath, ['elem.joint.',actual,'.Phi'],[1, 1], [1, Inf]);
            thetap = ncread(outpath, ['elem.joint.',actual,'.Omega'],[1, 1], [1, Inf]);
            % Flap angle and velocity
            beta = -ncread(outpath, ['elem.joint.',actual,'.Phi'],[2, 1], [1, Inf]);
            betap = ncread(outpath, ['elem.joint.',actual,'.Omega'],[2, 1], [1, Inf]);

            % LS fit for F and D
            ax = nexttile(tlo(fid-3));
            functionplot(ax,mbdyn_tab.t,xii,step,Color{ii},Marker{ii},'none',leg1)
            functionplot(ax,mbdyn_tab.t,delta_f,step,Color{nxt},Marker{nxt},'none',leg2)
            functionplot(ax,mbdyn_tab.t,delta_d,step,Color{ant},Marker{ant},'none',leg3)
            xlabel('t [s]','Interpreter','latex')
            ylabel('$\xi$, $\delta$ [rad]','Interpreter','latex')
            grid on
            legend('show','interpreter','latex')
            
            % LS fit for delta and theta
            ax = nexttile(tlo(fid-2));
            yyaxis left
            functionplot(ax,mbdyn_tab.t,theta,step,Color{ii},Marker{ii},'none',legt)
            ylabel('$\theta$[rad]')
            yyaxis right
            functionplot(ax,mbdyn_tab.t,delta_f,step,Color{nxt},Marker{nxt},'none',leg2)
            functionplot(ax,mbdyn_tab.t,delta_d,step,Color{ant},Marker{ant},'none',leg3)
            xlabel('t [s]','Interpreter','latex')
            ylabel('$\delta$ [rad]','Interpreter','latex')
            grid(ax,"on")
            legend(ax)

            % LS for delta and beta
            ax = nexttile(tlo(fid-1));
            yyaxis left
            functionplot(ax,mbdyn_tab.t,beta,step,Color{ii},Marker{ii},'none',legb)
            ylabel('$\beta$[rad]')
            yyaxis right
            functionplot(ax,mbdyn_tab.t,delta_f,step,Color{nxt},Marker{nxt},'none',leg2)
            functionplot(ax,mbdyn_tab.t,delta_d,step,Color{ant},Marker{ant},'none',leg3)
            xlabel('t [s]','Interpreter','latex')
            ylabel('$\delta$ [rad]','Interpreter','latex')
            grid(ax,"on")
            legend(ax)


            % Contribution to damper movement caused by LL
            % * delta D arm
            delta_xid = geom.Kxidelta*xii+Kxx*xii.^2;

            % Obtain residual contributions from theta and beta
            delta_d1 = delta_d-delta_xid;
            % Error
            delta_dcor = delta_d1-Ktt*theta.^2-Kbb*beta.^2-Ktb*theta.*beta;

            % delta F arm
            delta_f1 = delta_f-delta_xid;
            delta_fcor = delta_f1+Ktt*theta.^2+Kbb*beta.^2+Ktb*theta.*beta;

            % Known values (example data)
            x = theta;  % Replace with your x values
            y = beta;  % Replace with your y values
            f = delta_d1;  % Replace with your f values

            % Construct matrix A and vector b
            A = [x.^2', y.^2',(x.*y)'];
            b = f';

            % Solve for K1 and K2
            K = A \ b;

            % Extract K1 and K2
            K1 = K(1);
            K2 = K(2);
            K3 = K(3);

            % Display the results
            disp(['K1 = ', num2str(K1)]);
            disp(['K2 = ', num2str(K2)]);
            disp(['K3 = ', num2str(K3)]);

            % K(1) = 0.1;
            % K(2) = -0.2;
            % K(3) = 0.34;
            ax = nexttile(tlo(fid));
            colororder(ax,{'k','k'})
            yyaxis left
            plot(ax,mbdyn_tab.t,abs(delta_dcor),'ro-','DisplayName','Analytical'); hold on
            xlabel(ax,'$t$ [s]')
            ylabel(ax,'Abs. Error Analytical')
            yyaxis right
            plot(ax,mbdyn_tab.t,abs(delta_d1'-A*K),'bo-','DisplayName','LS fit')
            xlabel(ax,'$t$ [s]')
            ylabel(ax,'Abs. Error MBDyn')
            legend(ax)
            grid on
            xlim([mbdyn_tab.t(1) mbdyn_tab.t(end)])

            fid = fid+1;
            fhandle = figure(fid);
            fhandle.Name = ['Nb' int2str(Nb) damping{i} 'KINCOUP'];
            plot(mbdyn_tab.t,xii,'DisplayName','$\xi$'); hold on
            plot(mbdyn_tab.t,xii.^2,'DisplayName','$\xi^2$')
            plot(mbdyn_tab.t,beta.^2,'DisplayName','$\beta^2$')
            plot(mbdyn_tab.t,theta.^2,'DisplayName','$\theta^2$')
            plot(mbdyn_tab.t,theta.*beta,'DisplayName','$\theta\beta$')
            plot(mbdyn_tab.t,delta_d,'DisplayName','$\delta_D$')
            xlabel('$t$ [s]')
            ylabel('Angles')
            legend
            grid on
            xlim([mbdyn_tab.t(1) mbdyn_tab.t(end)])
            return
        end
        
        fhandle = figure(fid); 
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'damper'];
        tlo(fid) = tiledlayout('flow');
        for ii=1:Nb
            colororder({'k','k'})
            label = int2str(10000 + 1000*ii+30+1);
            xi_i = ncread(outpath,['elem.joint.' label '.E'],[3 1],[1 inf])*pi/180;
            % xi_i = xi_i(aux:end);
            xi_di = ncread(outpath,['elem.joint.' label '.Omega'],[3 1],[1 inf]);
            % xi_di = xi_di(aux:end);
            ax = nexttile(tlo(fid));
            yyaxis left
            functionplot(ax,mbdyn_tab.t,xi_i,step,Color{ii},Marker{ii},'none',['$\Delta\delta_' num2str(ii) '$'])
            ylabel('$\Delta\delta$ [rad]')
            yyaxis right
            functionplot(ax,mbdyn_tab.t,xi_di,step,Color{ii+1},Marker{ii},'none',['$\Delta\dot{\delta}_' num2str(ii) '$'])
            % plot(mbdyn_tab.t,xi_d,'DisplayName',['$\dot{\xi}_' num2str(ii) '$']); hold on
            ylabel('$\Delta\dot{\delta}$ [rad/s]')
            xlabel('t [s]')
            legend
            grid on

        end
        fid = fid+1;

        fhandle = figure(fid); fid = fid+1;
        fhandle.Name = ['Nb' int2str(Nb) damping{i} 'Md'];
        tlo(fid) = tiledlayout('flow');
        for ii=1:Nb
            defhinge = int2str(10000 + 1000*ii+30+1);
            md_i = ncread(outpath,['elem.joint.' defhinge '.m'],[3 1],[1 inf]);
            % xi_i = xi_i(aux:end);
            % xi_di = xi_di(aux:end);
            ax = nexttile(tlo(fid));
            functionplot(ax,mbdyn_tab.t,md_i,step,Color{ii},Marker{ii},'none',['$M_{d' num2str(ii) '}$'])
            ylabel('$M_d$ [Nm]')
            xlabel('t [s]')
            legend
            grid on
        end

    end


    
    %% Figures
    
    psi_nd = ncread(outpath,'time');
    % % #1, Lead-lag of single blade
    % % Figure = kk*10+1
    % % Tiled layout = kk*10+1
    fig(fid) = figure(fid);
    fig(fid).Name = [fignames 'llsb'];
    tlo(fid) = tiledlayout('flow'); lastfid = fid;
    fid = fid +1 ;

    psi = ncread(outpath, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    for ii = 1:Nb
        blade = int2str(10000+1000*ii + 30);
        xi(ii,:) = ncread(outpath, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);

        ax(1) = nexttile(tlo(lastfid));
        hold(ax(1),"on");
        functionplot(ax(1),psi_nd,-180/pi*xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        grid(ax(1),"on")
        xlabel(ax(1),'revolutions [-]');
        ylabel(ax(1),'$\xi$ [deg]')

    end

    % #9, Pith-flap lag
    joint_label = int2str(10000 + 1000*Nb + 30);
    label = joint_label;
    fig(fid)=figure(fid); 
    fig(fid).Name = [fignames 'PFL'];
    blade_txt = ['blade ', int2str(Nb)];
    fid = fid+1;
    % transform in deg, opposite sign for flap and lead-lag
    plot(psi_nd, 180/pi*ncread(outpath, ['elem.joint.',label,'.Phi'], [1, 1], [1, Inf]), ...
        psi_nd, -180/pi*ncread(outpath, ['elem.joint.',label,'.Phi'], [2, 1], [1, Inf]), ...
        psi_nd, -180/pi*ncread(outpath, ['elem.joint.',label,'.Phi'], [3, 1], [1, Inf]));
    xlabel('revolutions [-]');
    ylabel([blade_txt, ' angles, [deg]']);
    legend('$\theta$', '$\beta$', '$\xi$');
    grid on

    fid = fid+10;


    % clear t sol
    end
end
%%
% FolderName=strcat(pwd,'\images');
% FoldeName=char(FolderName);
% % FolderName = '\Images';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   saveas(FigHandle, fullfile(FolderName, [FigName, '.eps']),'epsc');    %<---- 'Brackets'
% end