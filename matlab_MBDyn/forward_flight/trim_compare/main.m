close all; clear all; clc
%% Color map and plot properties
%% Plotting parameters
setPlot
% cmap = colormap(lines(8));
% close all
% Color =cell(1,8);
% for ii=1:8
%    Color{ii}= cmap(ii,:);
% end
% Marker = {'o','+','x','s','d','^','*','.'};
% LineStyle = {'-' '--' ':' '-.' '-' '--' ':' '-.'};
% step = 128;
% plt.Color = Color; plt.Marker = Marker;
% plt.LineStyle= LineStyle;
% plt.step = step;
% 
% set(0,"defaultLineLineWidth",2,'DefaultAxesTickLabelInterpreter', 'latex',...
%     'DefaultTextInterpreter','latex','DefaultLegendInterpreter','latex');
fid = 1;
lastfid = fid;

%% Problem definition

ptype = 'inverse';
% Three prblem types can be defined
% direct: Input B1 and theta0 and solve trim
% inverse: Input Vinf and tau and solve trim
% autorot: Input Vinf and Q=0 and solve trim
%   * This problem is not a 14 equation system but a 15 equation system as
%   tau is also a parameter to be obtained as a solution.

[data,input,x0] = probdef(ptype);

% mu_vec = 0.:0.001:0.45;
% V = linspace(0.01,160,200);
% V = 50:5:80;
% V = [70 75 80];
V = [50:5:65 75 80] ;
% V = [65 75 80];

counter=1;
pid = 1;
for ii=1:length(V)
    input.V = V(ii);

    [x,fval,exitflag(ii),~]=fsolve(@(x)equations(x,ptype,input,data),x0);
    
    sol = extractdata(ptype,x,input);

    if ~strcmp(ptype,'autorot')
        sol.Q = data.rho*data.v_tip^2*data.A*data.R*data.sigma*(1/8*data.cd0*(1+sol.mu^2)+...
            1/2*data.cla*(sol.lambda*(1/3*sol.theta0-1/2*sol.lambda+sol.mu*(1/2*sol.a1-1/4*sol.B1))));
    end
    
    if exitflag(ii)>0
        V_inf(counter) = sol.V; V_z(counter) = sol.V*sin(sol.tau); 
        % alpha (counter) = sol.alfaD;
        % tau (counter)=sol.tau;
        sol.alfaH = sol.gamma+sol.tau;
        sol.muH = sol.V*cos(sol.alfaH)/(data.Omega*data.R);
        
        a1_H = sol.a1-sol.B1;
        sol.a1_H = a1_H*180/pi;
        sol.cT = sol.TD/(data.rho*(data.Omega*data.R)^2*data.A);
        sol.cQ = sol.Q/(data.rho*(data.Omega*data.R)^2*data.A*data.R);
        sol.TH = sol.TD*cos(a1_H)-sol.HD*sin(a1_H);
        sol.HH = sol.TD*sin(a1_H)+sol.HD*cos(a1_H); 
        sol.lambdai = sol.u/(data.Omega*data.R);
        % sol.B1 = sol.B1*180/pi;
    
        % if mod(counter,13)==0 && sol.mu<0.4
        trim{pid}=sol;
        pid=pid+1;
        % end
        counter=counter+1;
    end
end
clc
% % Assuming your cell array is named 'cellArray' containing structures with the same field names
% 
% % Initialize the combined structure
% trimstr = struct();
% 
% % Iterate through each structure in the cell array
% for i = 1:numel(trim)
%     % Get the current structure
%     currentStruct = trim{i};
% 
%     % Iterate through each field in the current structure
%     fields = fieldnames(currentStruct);
%     for j = 1:numel(fields)
%         % Get the field name
%         fieldName = fields{j};
% 
%         % Check if the field already exists in the combined structure
%         if isfield(trimstr, fieldName)
%             % Append the field value to the existing field
%             trimstr.(fieldName) = [trimstr.(fieldName), currentStruct.(fieldName)];
%         else
%             % Create a new field in the combined structure and assign its value
%             trimstr.(fieldName) = currentStruct.(fieldName);
%         end
%     end
% end
% % [results,fid] = rotortrim(data,trimstr,fid);
% fid = plotting(results,trimstr,fid,'','ff');
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%======================= GENERAL HELICOPTER parameters=================
mediumhelicopter
% Select number of blades

Nb = b;         % [-], number of blades

omega = Omega;     %[rad/s], rotor angular speed
nper = 20;      % [-], number of periods to run
% R = ;        % [m], rotor radius
% e = 1/25*R;     % [m], hinge position from hub
nuxi = sqrt(3/2*e/(R-e));

% Activate or deactivate additional blade-hub damper in blades
% 0 : inactive
% 1 : active
% C_b2h = b2hdamp*C_xi(gamma~0.8)
b2hdamp = 0.;
% Activate or deactivate blade degrees of freedom
% Since the degrees of freedom are given by means of a total joint in MBDyn
% 0 means that the degrre of freedom is not constrained and viceversa for
% 1.
% 0: allowed
% 1: clamped
flap = 0;
lag = 0;

% Execute or no MBDyn
% 1 : execute MBDyn
% else: do not execute
exec = 1;

% Geometrical consideration
opt = 1;

% Activate or deactivate swashplate dof
swp = 1;

% If the swashplate is active the pitch dof should be allowed again
if swp == 1
    pitch = 0;
else
    pitch = 1;
end

damp = {'std'};
%% MBDyn variables
geom = calculateKgeom(opt,e,damp,Nb);
%% MBDyn pathfiles

% path files for MBDyn
pref_in = '/home/ale/thesis/AERO_MBDyn/';
pref_out = ['/mnt/c/Users/aasal/Documents/Universidad/Master/Thesis/Codes/matlab_MBDyn/...' ...
    'forward_flight/trim_compare/'];
folder = 'aerodynamic_sw/';
damping = damp{1};
sf_out = 'v1/';

% if swp==1
%     filename = ['rotor_aero_' damping '_Nb' int2str(Nb) '_sw' '.mbd'];
%     outfile = [damping '_Nb' int2str(Nb) '_sw'];
% else
%     filename = ['rotor_aero_' damping '_Nb' int2str(Nb) '_nsw' '.mbd'];
%     outfile = [damping '_Nb' int2str(Nb) '_nsw'];
% end
if swp==1
    filename = ['aw169_' damping '_Nb' int2str(Nb) '_sw' '.mbd'];
    outfile = [damping '_Nb' int2str(Nb) '_sw'];
else
    filename = ['aw169_' damping '_Nb' int2str(Nb) '_nsw' '.mbd'];
    outfile = ['aw169' damping '_Nb' int2str(Nb) '_nsw'];
end

fignames = damping;
plt.fignames = fignames;
fid =fid;
% Input file path
in = [pref_in folder];
% Output file path seen from wsl pov
out = [pref_out sf_out];
% Output file path seen from the windows pov
outwd = sf_out;
% Execute MBDYN
outpath = [sf_out outfile '.nc'];
%% Different figures
lastfid = fid;
% #1
fh(fid) = figure(fid); 
fh(fid).Name = ['forces'];
tlo(fid) = tiledlayout('flow');

fid = fid+1;
% #2
fh = figure(fid); 
fh.Name = ['moments'];
tlo(fid) = tiledlayout('flow');

fid = fid+1;

% # 3
fh = figure(fid); 
fh.Name = ['vi'];
tlo(fid) = tiledlayout('flow');

fid = fid+1;

% #4
fh = figure(fid); 
fh.Name = ['mu'];
tlo(fid) = tiledlayout('flow');

fid = fid+1;

% #5
fh = figure(fid); 
fh.Name = ['inflow' ];
tlo(fid) = tiledlayout('flow');

fid = fid+1;

% #6
fig(fid) = figure(fid);
fig(fid).Name = [plt.fignames 'llallb'];
tlo(fid) = tiledlayout('flow');

fid = fid+1;
%%
for jj=1:length(trim)
    % Inputs for mbdyn
    V_INF= trim{jj}.V;
    ALTITUDE = 0.;
    INPUT_THETA_0=trim{jj}.theta0;
    INPUT_A_1_H = 0*pi/180; % radian
    INPUT_B_1_H=trim{jj}.B1;
    ALPHA_H=trim{jj}.alfaH;

variables =['"' 'const real R = ' num2str(R) ...
    '; const real M_BLADE =' num2str(Mb) ...
    '; const real S_BLADE_Z = ' num2str(Sb) ...
    '; const real J_BLADE_Z = ' num2str(Ib) ...
    '; const real X_BLADE_FLAP_OFFSET = ' num2str(e) ...
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
    mbdynsol{jj}=executembdyn(in,out,outwd,filename,outfile,variables,fid,plt,exec,V_INF,lastfid,tlo);
    fid=fid+10;
end
%%
for jj=1:length(trim)
trim_mat(:,jj) = cell2mat(struct2cell(trim{jj}));
mbdyn_mat(:,jj) = cell2mat(struct2cell(mbdynsol{jj}));
end
mbdyn_tab = array2table(mbdyn_mat','VariableNames',fieldnames(mbdynsol{1}));
trim_tab = array2table(trim_mat','VariableNames',fieldnames(trim{1}));
V_tab = array2table(trim_tab.V,'VariableNames',{'V'});
trim_tab = trim_tab(:,mbdyn_tab.Properties.VariableNames);
trim_tab.B1 = trim_tab.B1*180/pi;
% relerr_tab =table('VariableNames',fieldnames(mbdynsol{1}));
relerr_tab = abs((trim_tab-mbdyn_tab)./mbdyn_tab).*100;
mbdyn_tab = [V_tab mbdyn_tab];
trim_tab = [V_tab trim_tab];
% relerr_tab = [V_tab relerr_tab];
table2latex(relerr_tab,'relerrMBDyn')
table2latex(mbdyn_tab,'MBDynresults')
table2latex(trim_tab,'fsolveresults')
filename = 'errors.xlsx';
writetable(trim_tab,filename,'Sheet',1,'Range','A1')
writetable(mbdyn_tab,filename,'Sheet',1,'Range','A8')
writetable(relerr_tab,filename,'Sheet',1,'Range','A15')
%%
legend_ii ={'fsolve' 'MBDyn'};
n=0;
fig_handle=figure(fid+n); n = n+1;
fig_handle.Name ='TH';
plot(trim_tab.V,trim_tab.TH,'b-o'); hold on
plot(trim_tab.V,mbdyn_tab.TH,'r-o'); hold on
% errorbar(trim{jj}.V,trim{jj}.TH,abs(mbdynsol{jj}.T_H-trim{jj}.TH)); hold on
grid on
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')
ylabel('$T_H$ [N]','Interpreter','latex')
xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])
legend(legend_ii,'Interpreter','latex','Location','best')

fig_handle=figure(fid+n); n = n+1;
fig_handle.Name ='HH';
plot(trim_tab.V,trim_tab.HH,'b-o'); hold on
plot(trim_tab.V,mbdyn_tab.HH,'r-o'); hold on
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')
% errorbar(trim{jj}.V,trim{jj}.HH,abs(mbdynsol{jj}.H_H-trim{jj}.HH));hold on
grid on
xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])
ylabel('$H_H$ [N]','Interpreter','latex')
legend(legend_ii,'Interpreter','latex','Location','best')

% fig_handle=figure(fid+n); n = n+1;
% fig_handle.Name ='CH';
% plot(trim{jj}.V,trim_tab.Q,'b-o'); hold on
% plot(trim{jj}.V,mbdyn_tab.Q,'r-o'); hold on
%  % errorbar(trim{jj}.V,trim{jj}.Q,abs(mbdynsol{jj}.C_H-trim{jj}.Q));hold on
%  grid onn=0;
% xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])

fig_handle=figure(fid+n); n = n+1;
fig_handle.Name ='u';
plot(trim_tab.V,trim_tab.u,'b-o'); hold on
plot(trim_tab.V,mbdyn_tab.u,'r-o'); hold on
 % errorbar(trim{jj}.V,trim{jj}.u,abs(mbdynsol{jj}.u-trim{jj}.u));hold on
grid on
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')
ylabel('u [m/s]','Interpreter','latex')
xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])
legend(legend_ii,'Interpreter','latex','Location','best')

fig_handle=figure(fid+n); n = n+1;
fig_handle.Name ='t0';
plot(trim_tab.V,trim_tab.theta0,'b-o'); hold on
plot(trim_tab.V,mbdyn_tab.theta0,'r-o'); hold on
 % errorbar(trim{jj}.V,trim{jj}.theta0,abs(mbdynsol{jj}.theta_0-trim{jj}.theta0));hold on
grid on
ylabel('$\theta_0$ [deg]','Interpreter','latex')
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')
xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])
legend(legend_ii,'Interpreter','latex','Location','best')

fig_handle=figure(fid+n); n = n+1;
fig_handle.Name ='B1';
plot(trim_tab.V,trim_tab.B1,'b-o'); hold on
plot(trim_tab.V,mbdyn_tab.B1,'r-o'); hold on
 % errorbar(trim{jj}.V,trim{jj}.B1,abs(mbdynsol{jj}.B_1_H-trim{jj}.B1));hold on
 grid on
ylabel('$B_1$ [deg]','Interpreter','latex')
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')

xlim([round(trim_tab.V(1)) round(trim_tab.V(end))])
legend(legend_ii,'Interpreter','latex','Location','best')

fig_handle=figure(fid+n);
fig_handle.Name = 'relerr';
% for kk=1:size(relerr_tab,2)
bar(trim_tab.V,table2array(relerr_tab(:,1:7))'); hold on
% end
ylabel('Rel. error \%','Interpreter','latex')
ylim([0 100])
xlabel('$V_{\infty}$ [m/s]','Interpreter','latex')
xlim([round(trim_tab.V(1))-5 round(trim_tab.V(end))+5])
legend_ii = {'$T_H$','$H_H$','Q','u','$\mu_H$','$\theta_0$','$B_1$','$a_{1H}$'};
legend(legend_ii,'Interpreter','latex','Location','best')
grid on
return
%%
FolderName=strcat(pwd,'\images');
FoldeName=char(FolderName);
% FolderName = '\Images';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, fullfile(FolderName, [FigName, '.eps']),'epsc');    %<---- 'Brackets'
end