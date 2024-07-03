function [params, mbdynsol] = executembdyn(in,out,outwd,filename,outfile,fid,variables,exec)
%% Filpath definitions
% filepathin = '/home/ale/thesis/';
% filepathout = '/mnt/c/Users/aasal/Documents/Universidad/Master/Thesis/Codes/';
shfile = '/home/ale/mbdyn/mbdyn/mbdyn';

%% Exectue MBDyn

% Output file name
fn_base = outfile;

% Definition of enviroment variables for MBDyn
% variables =['"' 'const real V_INF = ' num2str(V_INF) ...
%     '; const real ALTITUDE = ' num2str(ALTITUDE) ...
%     '; const real ALPHA_H = ' num2str(ALPHA_H) ...
%     '; const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0) ...
%     '; const real INPUT_A_1_H = ' num2str(INPUT_A_1_H) ...
%     '; const real INPUT_B_1_H = ' num2str(INPUT_B_1_H) '"'];
system(['wsl MBDYNVARS=' variables]);
% system('echo %MBDYNVARS%')  

f_in = [in filename];
f_out = [out outfile];
disp('executing MBDyn...');
%%% WSL (might need to fix some paths...):
% [rc, errmsg] = system(['wsl MBDYNVARS="const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0) ' " ./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
% [rc, errmsg] = system(['wsl ./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
if exec == 1
    [rc, errmsg] = system(['wsl MBDYNVARS=' variables ' ' shfile ' ' f_in ' -o ' f_out]);
    if (rc ~= 0)
        error(errmsg);
    end
    disp('   ... done');
end
fn_outwd = [outwd fn_base];
% filename
fn = ['./' fn_outwd '.nc'];

disp(sprintf('reading from file ''%s''', fn));
%% Declaration of variables and numbering
GROUND = 0; 
AIRFRAME_X = 100; 
AIRFRAME_Y = 200;
HUB = 10000; 
 
BLADE_FLAP_OFFSET = 30; 
BLADE_CM = 40; 

A_NOD = 50;
B_NOD = 60;
C_NOD = 70;
D_NOD = 80;
F_NOD = 90;

% BLADE_TIP = 44;
% angular velocity
%%% Omega = 40; % rad/s
% we could get it from the motion of the hub, assuming it starts with nominal RPM
Omega = ncread(fn, ['node.struct.' num2str(HUB) '.Omega'], [3, 1], [1, 1]);
params.omega = Omega;
% compute azimuth vector from database
T = 0.;
time = ncread(fn,'time');
index = time>=0;
params.index = find(index==1,1);
time = time(index)-T;
dt = time(2)-time(1);
mbdynsol.t = time;
psi_nd = ncread(fn, ['elem.joint.' num2str(HUB) '.Phi'], [3, 1], [1, Inf]);
psi_nd = psi_nd(index);
% mbdynsol.psi_nd = psi_nd';
nT = 128;
revN = fix(length(psi_nd)/nT) - 1;
% extract number of blades from database
finfo = ncinfo(fn);
nb = 0;
while (1)
    nb = nb + 1;
    label = int2str(10000 + 1000*nb);
    n = 1;
    while (n <= size(finfo.Variables, 2))
        if (strcmp(finfo.Variables(n).Name, ['node.struct.',label,'.X']))
            break;
        end
        n = n + 1;
    end
    if (n > size(finfo.Variables, 2))
        break;
    end
end

nb = nb - 1;
params.nb = nb;
deltapsi = 2*pi/nb;
BLADE = zeros(1,nb);
for i = 1: nb
BLADE(i) = HUB+1000*i;
end
%% Post-processing
x = num2str(AIRFRAME_X);
y = num2str(AIRFRAME_Y);

hub_x = ncread(fn, ['node.struct.',x,'.X'],[1, 1], [1, Inf]);
hub_x = hub_x(index)';
hub_y = ncread(fn, ['node.struct.',y,'.X'],[2, 1], [1, Inf]);
hub_y = hub_y(index)';
% fighandle = figure(fid); fid = fid+1;
% fighandle.Name = 'Airframe';
% plot(time,hub_x); hold on
% plot(time,hub_y); hold on
% grid on

mbdynsol.hub_x= hub_x;
mbdynsol.hub_y= hub_y;
nxt_blade = circshift(BLADE,-1);
ant_blade = circshift(BLADE,1);
% Select blade number
for ii = 1:nb
ant = mod(ii - 2, nb) + 1;  % Calculate previous index
nxt = mod(ii, nb) + 1;  % Calculate next index
% node_label = int2str(BLADE(ii));
prev = int2str(BLADE(ant) + D_NOD);
actual = int2str(BLADE(ii) + BLADE_FLAP_OFFSET);
next = int2str(BLADE(nxt) + F_NOD);
psi = psi_nd+ii*deltapsi;
% psi = ncread(fn, ['node.struct.',node_label,'.E'],[3, 1], [1, Inf])*pi/180;
% rel_pos = pos()
xi = ncread(fn, ['elem.joint.',actual,'.Phi'],[3, 1], [1, Inf]);
% xi = aermech(3,:);
mbdynsol.(['xi_' num2str(ii)]) = xi(index)';
if contains(filename,'i2b')
    delta_xi = ncread(fn, ['elem.joint.',int2str(BLADE(ii) + BLADE_FLAP_OFFSET+1),'.E'],[3, 1], [1, Inf])-pi/2;
    delta_d = ncread(fn, ['elem.joint.',prev,'.Phi'],[3, 1], [1, Inf]);
    delta_f = ncread(fn, ['elem.joint.',next,'.Phi'],[3, 1], [1, Inf]);
    mbdynsol.(['Deltaxi_' num2str(ii)]) = delta_xi(index)';
    mbdynsol.(['delta_D' num2str(ant)]) = delta_d(index)';
    mbdynsol.(['delta_F' num2str(nxt)]) = delta_f(index)';
end
mbdynsol.(['psi_' num2str(ii)]) = psi';
% fighandle = figure(fid); %fid = fid+1;
% fighandle.Name = 'leadlag';
% plot(time,xi,'Marker','o','LineStyle','none','DisplayName',['$\xi_' num2str(ii) '$, MBDyn']); hold on
% legend('show','interpreter','latex')
% grid on

% joint_label = int2str(BLADE(ii)+BLADE_FLAP_OFFSET+1);

% % data = ncread()
% Y = fft(xi(nT*revN+[1:nT]));
% Y = Y/(nT/2);
% Y(1) = Y(1)/2;
% disp(sprintf('    a_0 = %g deg', Y(1)*180/pi));
% disp(sprintf('    a_1_H = %g deg', -real(Y(2))*180/pi));
% disp(sprintf('    b_1_H = %g deg', +imag(Y(2))*180/pi));
% mbdynsol.xi_1c = -real(Y(2))*180/pi;
% mbdynsol.xi_1s = +imag(Y(2))*180/pi;
end