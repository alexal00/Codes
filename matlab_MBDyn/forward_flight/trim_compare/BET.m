function [forces]=BET(data,lambdai_k)
% Extract data from the information structure data
% Control law
% theta = theta0+theta1*y/R+theta1s*sin(psi)+theta1c*cos(psi)
theta0=data.theta0; % Collective pitch, [rad]
thetaw=data.thetaw; % % Linear twist, [rad]
theta1c = data.theta1c; theta1s = data.theta1s; % Lateral and longitudinal pitch, [rad]

% Geometrical parameters of the rotor
R=data.R;       % Radius, [m]
y0=data.e;     % Cut-out/hinge position, [m]
c=data.c;       % Chord, [m]

y=data.y;       % Blade discretization
dy = y(2)-y(1);
y = y(y<=R*data.B); % Effective blade
r=y./R;         % N.D blade stations

cla=data.cla;   % Lift slope of airfoil
Omega=data.Omega; % Angular speed, [rad/s]
Ibeta = data.Ib;    % Inertia moment around hinge, [kgm^2]
Nb=data.b;     % Number of blades
mb=data.Mb/(data.R-data.e);     % mass per unit lenght of blades, [kg/m]
S = data.A;     % Rotor surface, [m^2]
% Flow porperties
rho=data.rho;   % Density, [kg/m^3]
a=data.a;       % Speed of sound, [m/s]
% RESULTS FROM THE TRIMMER
v = data.v;         % Freestream speed, [m/s]
alfaD = data.alfaD; % Angle of attack of the disk, [rad]
mu = data.mu ;      % N.D. airspeed, v/(OmegaR)
mux = data.mux;
muz = data.muz;
vi_k = lambdai_k*Omega*R;   % Induced velocity, [m/s]
if ~isscalar(vi_k)
    u = 1/R*(trapz(data.y,vi_k));       % Average induced velocity
    lambda=(v*sin(alfaD)+u)/(Omega*R);  % Inflow parameter
else
    u = vi_k;
    lambda=(v*sin(alfaD)+vi_k)/(Omega*R);
end
data.lambda = lambda;

lock=rho*c*cla*R^4/Ibeta; % lock no.
data.lock = lock;

% AEROMECHANICS EQUATIONS
a_ = aermec(data);
% Results from aeromechanics
% Beta = beta0+beta1c*cos(psi)+beta1s*sin(psi)
a0 = a_(1);     % Coning angle, [rad]
a1 = round(a_(2),3);     % Longitudinal flapping, [rad]
a2 = round(a_(3),3);     % Lateral flapping, [rad]
forces.a1 = a1;
forces.a2 = a2;

% 10 points & weights for Gauss integration along blade
% [nodes,weights]=lgwt(10);
nodes=[-0.9739065285171717,-0.8650633666889845,-0.6794095682990244,...
    -0.4333953941292472,-0.1488743389816312,0.1488743389816312,...
    0.4333953941292472,0.6794095682990244,0.8650633666889845,...
    0.9739065285171717];
weights=[0.0666713443086881,0.1494513491505806,0.2190863625159820,...
    0.2692667193099963,0.2955242247147529,0.2955242247147529,...
    0.2692667193099963,0.2190863625159820,0.1494513491505806,...
    0.0666713443086881];


% Values for numerical integration in psi with Gauss-Lobato quadrature
ninpsi=4;
dpsi=2.*pi/ninpsi;

% Initialise values for:
w = 0.;     % Power, [W]
q = 0.;     % Torque, [W]
Ft = 0.;    % Tangencial force in blades, [N]
T = 0.;     % Thrust, [T]
M_c = 0;    % Centrifugal moment at hinge, [Nm]
M_a = 0;    % Aerodynamic moment at hinge, [Nm]

for inpsi=1:ninpsi      % Loop for integration in psi
    psi1=(inpsi-1)*dpsi;
    psi2=inpsi*dpsi;
    % Integration in psi
    for i=1:length(nodes)
        psi=(psi1+psi2)/2.+nodes(i)*dpsi/2.;
        % Initialise values
        dTb = zeros(size(y));
        dFtb = zeros(size(y));
        dwr = zeros(size(y));
        dqr = zeros(size(y));
        dMb_a = zeros(size(y));
        dMb_c = zeros(size(y));
        % Integration in dr
        for k=1:length(y)
            rk=r(k);
            yk=y(k);
            % Induced velocity model
            if mu~=0
                % mux = mu*cos(alfaD);
                % Drees model for induced velocity
                chi=atan2(mux,lambda);
                kc=4/3*(1-cos(chi)-1.8*mu^2)/sin(chi);
                ks=-2*mu;
                vind=u*(1+rk*(kc*cos(psi)+ks*sin(psi)));
            else
                vind = vi_k(k);
            end
            % Harmonic flapping law
            Beta = a0+a1*cos(psi)+a2*sin(psi);
            Betadot = Omega*(-a1*sin(psi)+a2*cos(psi));
            % Pitching harmonic law
            theta = theta0+(thetaw*yk/R)+theta1c*cos(psi)+theta1s*sin(psi);

            % Airstream velocity components in TPP reference plane
            Vz=vind+v*sin(alfaD);   % Orthogonal to rotor disk
            Vy=v*cos(alfaD);        % Parallel to rotor disk         
            
            % Velocity components in the reference frame of the blade
            % UrP=Vy*cos(psi)-Vz*sin(Beta); % Radial velocity
            UrP = Vy*cos(psi)-Vz*Beta;
            UtP=Omega*yk+Vy*sin(psi); % Parallel to blade chord
            % UpP=Vz*cos(Beta)+Betadot*(yk-y0)+Vy*cos(psi)*sin(Beta); % Perpendicular velocity
            UpP=Vz+Betadot*(yk-y0)+Vy*cos(psi)*Beta; % Perpendicular velocity

            phi=atan2(UpP,UtP);
            GAMMA=atan2(UrP,UtP);
            
            alfa=theta-phi; % Local blade AoA
            U=sqrt(UtP^2+UpP^2);    % Total velocity, \approx UT
            
            Mach=U/a;
            [cl,cd,~]=cpcrcm(alfa,Mach);

            % Dynamic stall correction (Boeing Vertol method)
            % adot1 = (-theta1c*sin(psi)+theta1s*cos(psi))*Omega;
            % adot2 = 0;
            adot=(-theta1c*sin(psi)+theta1s*cos(psi))*Omega;
            if adot~=0
                if(Mach<.6)
                    Mm=Mach;
                elseif(Mach<=.3)
                    Mm=.3;
                end
                dalfa=61.5*log(.6/Mm)*sqrt(abs(adot)*c/2/abs(UtP));
                dalfa=dalfa*pi/180.*adot/abs(adot);
                [clst,~,~]=cpcrcm(alfa+dalfa,Mach);
                cl=clst*alfa/(alfa+dalfa);
            end

            % Radial velocity orrection (sweep angle)
            if ((cl/cos(GAMMA))<(cla*alfa))
                cl=cl/cos(GAMMA);
            end
            
            % Aerodynamic forces
            dL=.5*rho*U^2*cl*c;                 % Lift
            dD=.5*rho*(UrP^2+UtP^2+UpP^2)*cd*c; % Drag
            
            % Forces in blade reference frame
            dTb(k)=dL*cos(phi)-dD*sin(phi);     % Thrust
            dFtb(k)=dL*sin(phi)+dD*cos(phi)/cos(GAMMA); % Tangencial force
            % Power
            dwr(k)=dFtb(k)*UtP;
            dqr(k)=dFtb(k)*yk;
            % Moments
            dMb_a(k) = dTb(k)*(yk-y0);
            dMb_c(k) = mb*yk*(yk-y0)*Omega^2;
        end
        % Numerical integration in r
        % dT = trapz(y,dTb);
        dT = sum(dTb*dy);
        dFt = sum(dy*dFtb);
        dW = sum(dy*dwr);
        % dQ = trapz(y,dqr);
        dQ = sum(dqr*dy);
        dM_a = sum(dy*dMb_a);
        dM_c = sum(dy*dMb_c);
        % Numerical integration in psi
        M_a=M_a+dM_a*weights(i)*dpsi/2.;
        M_c=M_c+dM_c*weights(i)*dpsi/2.;
        w = w+dW*weights(i)*dpsi/2.;
        q = q+dQ*weights(i)*dpsi/2.;
        Ft = Ft+dFt*weights(i)*dpsi/2;
        T=T+dT*weights(i)*dpsi/2.;
    end      
end
% Multiply by the number of blades and average
Pbet=w*Nb/(2*pi);
% Qbet = Pbet / Omega;
Qbet=q*Nb/(2*pi);

Tbet=T*Nb/(2*pi);

M_a = Nb/(2*pi)*M_a;
M_c = Nb/(2*pi)*M_c;

% N.D. thurst and torque coefficients
CTbet = Tbet/(rho*(Omega*R)^2*S);
% % CQmom = Qmom/(rho*(Omega*R)^2*S*R); 
CPbet = Qbet/(rho*(Omega*R)^2*R*S);

Beta02 = M_a/M_c*180/pi;

% Store results in structure
forces.Beta0=a0*180/pi;
forces.Beta02=Beta02;
forces.T=Tbet; forces.CT=CTbet;
forces.P=Pbet; forces.CQ=CPbet;
forces.Q=Qbet;
end