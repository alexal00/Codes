function [coefs] = aermec(data) 
% Ode integration parameters
tspan = [0 10*2*pi];
beta0 = [0 0];
% ODE resolution
[psi , beta] =ode45(@(psi,x) odefun(psi,x,data),tspan,beta0);
% figure
% plot(psi,beta(:,1)); hold on
% plot(psi,beta(:,2));
nT = fix(length(psi)/10);
revN = fix(length(psi)/nT) - 1;
pt = nT*revN+[1:nT];
L = length(pt);
% Compute the FFT
Beta = fft(beta(pt,1));
% a0 = mean(beta(pt,1));
[a,b] = fourierFFT(3,Beta,L,1,1);

coefs = [a(1),a(2),b(2)];
end
function dydt = odefun(psi,x,data)
% EXAMPLE DATA
% c = 0.2;
% cla = 2*pi;
% theta0 = 3*pi/180;
% theta1 = 0;
% theta1c = 0*pi/180;
% theta1s = -2*pi/180;
% Omega =30;
% R = 5;
% rho = 1.225;
% lambda = 0.1;
% mu = 0.1;
% Ib = 100;
% e =0.1*R;
% nu_b = 1.22;
% y = linspace(0.1*R,R,20);
% Extract data from structure data
c = data.c; cla=data.cla; 
Omega = data.Omega; R = data.R;
rho = data.rho; lambda = data.lambda; mu = data.mu;
mux = data.mux;
a = data.a;
Ib = data.Ib; %Sb = data.Sb;
mb = data.Mb/(data.R-data.e);
e = data.e;
nu_b = data.nu_b;
y = data.y(data.y<=data.R*data.B);
% Harmonic pitching law in [rad]
theta0=data.theta0; theta1 = data.thetaw; theta1s = data.theta1s;
theta1c = data.theta1c; 

% Initialise values
mbeta=zeros(size(y));
% Integration in r
for i=1:length(y)
    yi = y(i);
    UT = Omega*yi+mux*Omega*R*sin(psi);
    Up = lambda*Omega*R+(yi-e)*Omega*x(2)+mux*Omega*R*x(1)*cos(psi);
    Ur = mux*Omega*R*cos(psi);
    U = sqrt(UT^2+Up^2);
    theta = theta0+theta1*yi/R+theta1c*cos(psi)+theta1s*sin(psi);
    phi = atan2(Up,UT);
    alfa = theta-phi;
    % Aerodynamic forces
    Mach=U/a;
    [cl,cd,~]=cpcrcm(alfa,Mach);
    dL=.5*rho*U^2*cl*c;                 % Lift
    dD=.5*rho*(Ur^2+UT^2+Up^2)*cd*c; % Drag
    L = 0.5*rho*U^2*c*cla*alfa;
            
    % Forces in blade reference frame
    dTb2=dL*cos(phi)-dD*sin(phi);     % Thrust
    dTb=L;     % Thrust
    mbeta(i) = dTb*(yi-e);
    mbeta2(i) = dTb2*(yi-e);
    Ib_i(i) = mb*yi^2;
end
% Ib = trapz(y,Ib_i);
% Aerodynamic moment at the flapping hinge
Mb = 1/(Ib*Omega^2)*trapz(y,mbeta);
% Aerodynamic moment at the flapping hinge
Mb2 = 1/(Ib*Omega^2)*trapz(y,mbeta2);
% ODE system
dydt = zeros(2,1);
dydt(1) = x(2);
dydt(2) = -nu_b^2*x(1)+Mb;
end