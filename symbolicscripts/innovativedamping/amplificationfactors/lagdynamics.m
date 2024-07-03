% Author: Alejandro Alvaro, 2023-2024
% evaluate frequency content of different arrangements
% under ideal excitation conditions and visualize ouput.
%
% Solutions are obtained after integrating the corresponding ODE system of
% equations related to the lead-lag motion of the blades, being excited by
% a periodic forcing term.
close all; clear all; clc
%% Plotting parameters
setPlot
%% Data for the ODE
omega = 40;             % Rotor angular speed, [rad/s]
T = 2*pi/omega;         % Revolution period, [s]

% Load data from the helicopter
data = hammond1974data;
T_xi = T/data.nu_xi;    % Lead-lag period, [s]

nT = 20;               % Number of lead-lag periods to run simulation
steps = 256;            % Steps per period
dtxi = T_xi/steps;      % time-step per period


dt = T/steps;
tspan = 0:dt:T*nT-dt;   % Time vector for ODE

Nb = 4;                 % Number of blades

% ODE initial conditions, [xi_i xidot_i]
y0 = [0*ones(1,Nb) zeros(1,Nb)];


% Select damping type
% 0 : Blade to hub
% 1 : Interblade
% 2 : Inter-2-blade
mode = [0 2];
name = {'std' 'i2b'};
% mode = 2;
% Maximum exciting harmonic in the force
n = 7;

lastfid = 1;

% #5
fig(100) = figure(100);
fig(100).Name = 'PSDDeltaxidotcomp';
for ii = 1:Nb
    ax(100+ii) = subplot(2,2,ii);
    hold on
end
% ax(100) = axes(fig(100));
% tlo(100) = tiledlayout('flow');

for i=1:length(mode)

    disp(['Damping, mode= ' num2str(mode(i)) ])
    disp(['Excitation harmonic, n= ' num2str(n)])
    % Solution of the ODE system
    if exist(['ld_mode' num2str(mode(i)) 'n' num2str(n) '.mat'],"file")
        load(['ld_mode' num2str(mode(i)) 'n' num2str(n) '.mat'],"sol","t")
    else
        [t,sol] = ode45(@(t,y) odefun2(t,y,omega,mode(i),Nb,n),tspan,y0);
        save(['ld_mode' num2str(mode(i)) 'n' num2str(n) '.mat'],"sol","t")
    end

    %% Plotting
    revN = fix(length(t)/steps)-1;

    % N.D time in revolutions
    revs = t/T;

    % revsxi = t/T_xi;

    % Select blade for plotting
    blades = [1 2 3 4];

    % Plot the lead-lag evolution of blade jj
    % #1
    fig(fid) = figure(fid); 
    fig(fid).Name = ['ximode' num2str(mode(i))];
    num_functions = length(blades);
    tlo(fid) = tiledlayout(num_functions, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    fid=fid+1;
    % #2
    fig(fid) = figure(fid); 
    fig(fid).Name = ['xidotmode' num2str(mode(i))];
    tlo(fid) = tiledlayout('flow');
    fid=fid+1;
    % #3
    fig(fid) = figure(fid);
    fig(fid).Name = ['Deltaxidotmode' num2str(mode(i))];
    tlo(fid) = tiledlayout('flow');
    fid=fid+1;

    % #4
    fig(fid) = figure(fid);
    fig(fid).Name = ['PSDDeltaxidotmode' num2str(mode(i))];
    tlo(fid) = tiledlayout('flow');
    % fid=fid+1;
    for jj = 1:length(blades)
        % Load the solution information
        xi = sol(:,jj);         % Lead-lag angle of blade jj
        xi_d = sol(:,Nb+jj);    % Lead-lag angular speed of blade jj
        
        ax(1) = nexttile(tlo(lastfid));
        % Function to calculate RMS
        xi_m = movmean(xi,steps);
        calculate_rms = @(funct) sqrt(mean(funct.^2));
        
        % Calculate RMS for each time step
        rms_values = arrayfun(@(n) calculate_rms(xi_m(1:n)), 1:length(xi));
        
        % Plot RMS values
        plot(ax(1),revs, -180/pi*xi,'DisplayName','$\xi$'); hold on
        plot(ax(1),revs,-180/pi*xi_m,'DisplayName','$\bar{\xi}$'); hold on
        xlabel(ax(1),'Rotor revolutions [-]');
        ylabel(ax(1),['$\xi_' num2str(blades(jj)) '$, mode' num2str(mode(i))]);
        legend
        grid on
        ax(2)=nexttile(tlo(lastfid));
        plot(ax(2),revs,rms_values); hold on
        ylabel(ax(2),'RMS');
        % legend
        grid on

        % Plot the lead-lag angular speed evolution of blade jj
        ax(3) = nexttile(tlo(1+lastfid));
        plot(ax(3),revs, -180/pi*xi_d); hold on
        xlabel('revolutions [-]');
        ylabel(['$\dot{\xi}_' num2str(blades(jj)) '$ [deg/s]']);
        legend
        grid on

        ax(4) = nexttile(tlo(2+lastfid));
        if mode(i)==0
            nxt = 0;
        elseif mode(i)==1
            nxt = mod(jj, Nb) + 1;  % Calculate next blade
        elseif mode(i) ==2
            nxt = mod(jj+1, Nb) + 1;  % Calculate next blade
        end
        if nxt ==0
            deltaxidot = xi_d;
        else
            xi_dnxt = sol(:,Nb+nxt);
            deltaxidot = xi_dnxt-xi_d;
        end
        plot(ax(4),revs,-180/pi*deltaxidot)
        xlabel(ax(4),'revolutions [-]');
        ylabel(ax(4),['$\Delta\dot{\xi}_' num2str(blades(jj)) '$ [deg/s]']);
        grid on

        % PSD of xi
        % fig(4) = figure(fid); fid=fid+1;
        % tlo(4) = tiledlayout('flow');
        % PSD of Deltaxi_dot
        %Compute fft and the evolution of the harmonics with time
        span = revN:revN;
        for kk=span
            if kk==span(1)
                time = t(1:steps);
            end

            % xi_aux= xi(steps*kk+[1:steps],:);
            xid_aux = deltaxidot(steps*kk+[1:steps],:);

            % FFT of the lead-lag angle
            % [P1xi,f]=computeFFT(xi_aux,time);

            % FFT of lead-lag angular speed
            [P1xid,f]=computeFFT(xid_aux,time);

            % Plot one-sided espectrum of lead-lag angular speed
            ax(5) = nexttile(tlo(3+lastfid));
            stem(ax(5),f,P1xid,'filled'); hold on
            xlabel(ax(5),'$f$ [Hz]')
            ylabel(ax(5),['$PSD(\Delta\dot{\xi}_' num2str(blades(jj)) ')$'])
            xlim(ax(5),[0 ceil(10*omega/(2*pi))])
            title(ax(5),['Revolution=' num2str(kk)])
            % legend
            grid on

            % Plot one-sided espectrum of lead-lag angular speed
            stem(ax(100+jj),f,P1xid,'filled',Color=Color{i},DisplayName=name{i}); hold on
            xlabel(ax(100+jj),'$f$ [Hz]')
            ylabel(ax(100+jj),['$PSD(\Delta\dot{\xi}_' num2str(blades(jj)) ')$'])
            xlim(ax(100+jj),[0 ceil(10*omega/(2*pi))])
            title(ax(100+jj),['Revolution=' num2str(kk)])
            legend(ax(100+jj))
            grid on

        end
    end
    % lastfid = fid;
    fid = i*10;
    lastfid = fid;
end

%% Auxiliary functions
function dydt = odefun2(time,y,omega,mode,blades,harmonic)
% Helicopter data
data = hammond1974data;
data.n_b = blades;
data.omega = omega;
% syms t real
subs_var = [data.S_b data.n_b data.m_b data.omega data.I_b ...
            data.c_xi data.k_xi data.e];
[M_R, C_R, K_R, sym_var] = lag_eq(mode,blades);
M_R = double(subs(M_R,sym_var, subs_var));
% M_R = matlabFunction(M_R);
% M_R = M_R(time);
C_R = double(subs(C_R,sym_var, subs_var));
% C_R = matlabFunction(C_R);
% C_R = C_R(time);
K_R = double(subs(K_R, sym_var, subs_var));
% K_R = matlabFunction(K_R);
% K_R = K_R(time);

dof = data.n_b;
dpsi = 2*pi/data.n_b;
xi_c = zeros(dof,1);
xi_s = zeros(dof,1);
f_xi = zeros(dof,1);
xi_0 = 0.01;
for ii=1:dof
    for jj=1:harmonic
    xi_c(ii) = xi_c(ii)+ xi_0*cos(jj*((omega*time)+(ii-1)*dpsi));
    xi_s(ii) = xi_s(ii)+ xi_0*sin(jj*((omega*time)+(ii-1)*dpsi));
    end
end
for jj=1:harmonic
    f_xi=f_xi+(-(omega*jj)^2*M_R+K_R)*xi_c-(omega*jj)*C_R*xi_s;
end

% f_xi=zeros(dof,1);
% for ii=1:dof
%     f_xi(ii)=xi_0*cos(omega*time+(ii-1)*dpsi);
% end

B = [zeros(dof);M_R\eye(dof)];
A = [zeros(dof) eye(dof);...
     -M_R\K_R -M_R\C_R];
dydt = A*y+B*f_xi;
aux = 1;
end