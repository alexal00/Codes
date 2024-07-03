function mbdynsol=executembdyn(in,out,outwd,filename,outfile,variables,fid,plt,exec,Vinf,lastfid,tlo)

%% File paths for MBDyn
% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn';

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

f_in = [in filename];
f_out = [out outfile];

disp('executing MBDyn...');
%%% WSL (might need to fix some paths...):
% [rc, errmsg] = system(['wsl MBDYNVARS="const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0) ' " ./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
% [rc, errmsg] = system(['wsl ./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
if exec ==1
    [rc, errmsg] = system(['wsl  MBDYNVARS=' variables ' ' shfile ' ' f_in ' -o ' f_out]);
    if (rc ~= 0)
        error(errmsg);
    end
    disp('   ... done');
end
fn_outwd = [outwd fn_base];
% filename
fn = [fn_outwd '.nc'];
% disp(sprintf('reading from file ''%s''', fn));

% angular velocity
%%% Omega = 40; % rad/s
% we could get it from the motion of the hub, assuming it starts with nominal RPM
Omega = ncread(fn, 'node.struct.10000.Omega', [3, 1], [1, 1]);

% compute azimuth vector from database
psi_nd = ncread(fn, 'time')*Omega/(2*pi);

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
deltapsi = 2*pi/nb;
blade_label = int2str(10000 + 1000*nb);
joint_label = int2str(10000 + 1000*nb + 30);
inflow_label = '99';

% disp(sprintf('  found %d blades; using blade %d, label=''%s''', nb, nb, blade_label));

% fid = 10;

label = inflow_label;
ax(1) = nexttile(tlo(lastfid));
plot(ax(1),psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.f']));
xlabel(ax(1),'revolutions, n.d.');
ylabel(ax(1),'rotor force, N');
legend(ax(1),'Longitudinal', 'Lateral', 'Thrust');
title(ax(1),['V = ' int2str(Vinf)])


nT = 256;
revN = fix(length(psi_nd)/nT) - 1;

F_H = ncread(fn, ['elem.inducedvelocity.',label,'.f']);
T_H = F_H(3, :);

Y = fft(T_H(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    T_H(0) = %g N', Y(1)));
mbdynsol.TH = Y(1);
disp(sprintf('    T_H(%d) = %g %+gi N', nb, real(Y(nb+1)), -imag(Y(nb+1))));

H_H = F_H(1, :);
Y = fft(H_H(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    H_H(0) = %g N', Y(1)));
mbdynsol.HH = Y(1);
disp(sprintf('    H_H(%d) = %g %+gi N', nb, real(Y(nb+1)), -imag(Y(nb+1))));

% fh = figure(fid); fid = fid+1;
% fh.Name = ['momentsV' int2str(Vinf)];
ax(1) = nexttile(tlo(lastfid+1));
plot(ax(1),psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.m']));
xlabel(ax(1),'revolutions, n.d.');
ylabel(ax(1),'rotor moments, Nm');
legend(ax(1),'Roll', 'Pitch', 'Torque');
title(ax(1),['V = ' int2str(Vinf)])

M_H = ncread(fn, ['elem.inducedvelocity.',label,'.m']);
C_H = M_H(3, :);
Y = fft(C_H(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    C_H(0) = %g Nm', Y(1)));
mbdynsol.Q = -Y(1);
disp(sprintf('    C_H(%d) = %g %+gi Nm', nb, real(Y(nb+1)), -imag(Y(nb+1))));

% fh = figure(fid); fid = fid+1;
% fh.Name = ['viV' int2str(Vinf)];
ax(1) = nexttile(tlo(lastfid+2));
plot(ax(1),psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.UMean']));
xlabel(ax(1),'revolutions, n.d.');
ylabel(ax(1),'mean induced velocity, m/s');
legend(ax(1),'u');
title(ax(1),['V = ' int2str(Vinf)])
ylim(ax(1),[0 1.1*max(ncread(fn, ['elem.inducedvelocity.',label,'.UMean']))])

UMean = ncread(fn, ['elem.inducedvelocity.',label,'.UMean']);
Y = fft(UMean(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    UMean(0) = %g m/s', Y(1)));
mbdynsol.u = Y(1);
% disp(sprintf('    UMean(%d) = %g %+gi m/s', nb, real(Y(nb+1)), -imag(Y(nb+1))));

mu = ncread(fn, ['elem.inducedvelocity.',label,'.Mu']);
Y = fft(mu(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
mbdynsol.muH=Y(1);
% fh = figure(fid); fid = fid+1;
% fh.Name = ['muV' int2str(Vinf)];
ax(1) = nexttile(tlo(lastfid+3));
plot(ax(1),psi_nd, mu);
xlabel(ax(1),'revolutions, n.d.');
ylabel(ax(1),'$\mu$, n.d.');
legend(ax(1),'advance ratio');
title(ax(1),['V = ' int2str(Vinf)])
ylim(ax(1),[0 1.1*max(ncread(fn, ['elem.inducedvelocity.',label,'.Mu']))])


lambda = ncread(fn, ['elem.inducedvelocity.',label,'.Lambda']);
Y = fft(lambda(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
% mbdynsol.lambda=Y(1);
% fh = figure(fid); fid = fid+1;
% fh.Name = ['inflowV' int2str(Vinf)];
ax(1) = nexttile(tlo(lastfid+4));
plot(ax(1),psi_nd, lambda);
xlabel(ax(1),'revolutions, n.d.');
ylabel(ax(1),'$\lambda$, n.d.');
legend(ax(1),'inflow ratio');
title(ax(1),['V = ' int2str(Vinf)])
ylim(ax(1),[0 1.1*max(lambda)])


label = joint_label;

psi_nd = ncread(fn, ['elem.joint.',label,'.Phi']);

theta = psi_nd(1, :);
Y = fft(theta(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    theta_0 = %g deg', Y(1)*180/pi));
disp(sprintf('    A_1_H = %g deg', -real(Y(2))*180/pi));
disp(sprintf('    B_1_H = %g deg', +imag(Y(2))*180/pi));
mbdynsol.theta0 = Y(1);
mbdynsol.B1 = +imag(Y(2))*180/pi;
beta = -psi_nd(2, :);
Y = fft(beta(nT*revN+[1:nT]));
Y = Y/(nT/2);
Y(1) = Y(1)/2;
disp(sprintf('    a_0 = %g deg', Y(1)*180/pi));
disp(sprintf('    a_1_H = %g deg', -real(Y(2))*180/pi));
disp(sprintf('    b_1_H = %g deg', +imag(Y(2))*180/pi));
mbdynsol.a1_H = -real(Y(2))*180/pi;
% mbdynsol.b1_H = +imag(Y(2))*180/pi;
% #1
% fig(fid)= figure(fid);
% fig(fid).Name = [plt.fignames 'llsb'];
% tlo(fid) = tiledlayout('flow');
% lastfid = fid;
% fid = fid +1 ;
% #2
% fig(fid) = figure(fid);
% fig(fid).Name = [plt.fignames 'llallb'];
ax(2) = nexttile(tlo(lastfid+5));
xlabel(ax(2),'revolutions [-]');
ylabel(ax(2),'$\xi$ [deg]')
hold on;
fid = fid +1 ;

psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
% compute azimuth vector from database
psi_nd = ncread(fn, 'time')*Omega/(2*pi);

for ii = 1:nb
    % ax(1) = nexttile(tlo(lastfid));
    % hold(ax(1),"on");
    psi_b(ii,:) = psi+ii*deltapsi;
    blade = int2str(10000+1000*ii + 30);
    xi(ii,:) = ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
    xid(ii,:) = ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);

    % functionplot(ax(1),psi_nd,-180/pi*xi(ii,:),plt.step,plt.Color{ii},plt.Marker{ii},plt.LineStyle{ii},['$\xi_' int2str(ii) '$'])
    % grid(ax(1),"on")
    % xlabel(ax(1),'revolutions [-]');
    % ylabel(ax(1),'$\xi$ [deg]')
    
    functionplot(ax(2),psi_nd,-180/pi*xi(ii,:),plt.step,plt.Color{ii},plt.Marker{ii},plt.LineStyle{ii},['$\xi_' int2str(ii) '$'])
    % plot(psi_nd,-180/pi*xi(ii,:),'DisplayName',['$\xi_' int2str(ii) '$']); 
end

end