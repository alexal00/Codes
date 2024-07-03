function rotorstate(data,index,name,tlo)
%
muz = data.muz;
mux = data.mux;
Omega = data.Omega;
R = data.R;
y0 = data.y0;
beta0 = data.Beta0;
beta1c = data.Beta1c;

Vz = [0; 0; -muz*Omega*R];
Vx = [0; 0; mux*Omega*R];
V = Vx+Vz+abs(Vz);

h = 0.5;
o = [0;0;0];
A = [0;0;h];
rco = [y0;0;0];
r1 = -rco+A;
r2 = rco+A;
beta = beta0+beta1c;
TAP =[cos(beta) 0 -sin(beta);...
        0 1 0;...
        sin(beta) 0 cos(beta)];
rt2=r2+TAP*[(R-y0);0;0];

beta = beta0-beta1c;
TAP =[cos(beta) 0 sin(beta);...
        0 1 0;...
        -sin(beta) 0 cos(beta)];
rt1=r1+TAP*[-(R-y0);0;0];
mast = [o,A];
disk = [r1,r2];
blade1 = [r1,rt1];
blade2 = [r2,rt2];
tpp = [rt1,rt2];

ax = nexttile(tlo);
plot(ax,mast(1,:),mast(3,:),'k'); hold on
plot(ax,disk(1,:),disk(3,:),'ko-'); hold on
plot(ax,blade1(1,:),blade1(3,:),'k');hold on
plot(ax,blade2(1,:),blade2(3,:),'k'); hold on
plot(ax,tpp(1,:),tpp(3,:),'r');hold on
title(ax,['Experiment ' num2str(index)])
xlim(ax,[-R R]);
ylim(ax,[0 0.7]);
end