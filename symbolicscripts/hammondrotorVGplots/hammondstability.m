close all; clear variables; clc;
%%
setPlot
%% Load helicopter data for Hammond rotor
% Select modes to plot the vg plot
% 0: std blade-to-hub
% 1: ib
% 2: i2b
mode = [0 1 2];
% Select whether or not to apply the generalized Deutsch criterion
deutsch = [1 2*(1-cos(2*pi/4)) 2*(1-cos(2*2*pi/4))];
% deutsch = [1 1 1];

% idx for internal loop
idx = 1;

%% Obtain coleman plot
for m = mode
    % Load data from Hammond rotor
    data = hammond1974data;

    % Safety factor for chosen lead lag damping
    data.SF = 1.;
    % Assumed aerodynamic parameters
    data.omega = 40;
    data.R = 5.5;
    % Correct the lead-lag damping to obtain same GR stability margins
    data.c_xi = data.c_xi./deutsch(idx);
    % Execute the code
    ae71(data,m);
    % Update the index of the loop
    idx = idx+1;
end

%%
return
savefigures