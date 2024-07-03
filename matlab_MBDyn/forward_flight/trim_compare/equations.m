function F= equations(x,ptype,input,data)
% data extracted from structure
% flow
rho = data.rho;
% construction
%general data
W = data.W ; sigma = data.sigma;
b = data.b ; Omega = data.Omega; R = data.R;
Sb = data.Sb; Kb = data.Kb; e = data.e;
A = data.A; cla = data.cla; cd0 = data.cd0;
%fuselage
f = data.f; cmf = data.cmf; Sbar = data.Sbar; lbar = data.lbar;
%tail
clt = data.clt; St=data.St;
%construction
xcg = data.xcg; l = data.l ; h=data.h;

if strcmp(ptype,'inverse')
    % V and tau are fixed
    V = input.V; tau = input.tau;
TD=x(1); HD=x(2);Rf=x(3); Mf=x(4); Pc=x(5);
a1=x(6); gamma=x(7); mu=x(8); lambda=x(9);
alfaD=x(10); u=x(11); v1=x(12);
theta0=x(13); B1=x(14);
% TD=x(1); HD=x(2);
% a1=x(3); gamma=x(4); mu=x(5); lambda=x(6);
% alfaD=x(7); u=x(8); v1=x(9);
% theta0=x(10); B1=x(11);
elseif strcmp(ptype,'direct')
    % theta0 and B1 are fixed
    theta0 = input.theta0; B1 = input.B1;
TD=x(1); HD=x(2); Rf=x(3); Mf=x(4); Pc=x(5);
a1=x(6); gamma=x(7); mu=x(8); lambda=x(9);
alfaD=x(10); u=x(11); v1=x(12);
V=x(13); tau=x(14);
else
    V = input.V; Q = input.Q;
TD=x(1); HD=x(2);Rf=x(3); Mf=x(4); Pc=x(5);
a1=x(6); gamma=x(7); mu=x(8); lambda=x(9);
alfaD=x(10); u=x(11); v1=x(12);
theta0=x(13); B1=x(14); tau = x(15);
end
% tip-speed
v_tip = Omega*R;
Kh=b/2*(e*Sb*Omega^2+Kb);
% Definition
% Vertical force of rotor
F(1)=TD-rho*sigma*A*v_tip^2*cla*1/2*(theta0*(1/3+mu^2/2)-lambda/2-mu*a1*0.5);
% horizontal force of rotor
F(2)=HD-rho*sigma*A*v_tip^2*(mu*cd0*0.25+cla*0.25*(mu*lambda*theta0-lambda*a1*0.5));
% Fuselage resistance
F(3)=Rf-0.5*rho*V^2*f;
% Fuselage moment
F(4)=Mf-0.5*rho*V^2*cmf*Sbar*lbar;
% Vertical tail
F(5)=Pc+0.5*rho*V^2*clt*St*(gamma+tau);
% Flapping longitudinal component
F(6)=a1-1/(1-0.5*mu^2)*(2*mu*(4/3*theta0-lambda));
% Advance ratio
F(7)=mu-V*cos(alfaD)/v_tip;
% inflow velocity
F(8)=lambda-(V*sin(alfaD)+u)/v_tip;
% induced velocity
F(9)=u-TD/(2*rho*A*v1);
% velocity in the plane of the disk
F(10)=v1-v_tip*sqrt(mu^2+lambda^2);
% angle of attack of the disk
F(11)=-alfaD+tau+gamma-a1+B1;
% Definition of the equilibirum equations
F(12)=TD*cos(a1-B1)-HD*sin(a1-B1)-W*cos(gamma)-Rf*sin(gamma+tau)+Pc*cos(gamma+tau);
% F(9)=TD*cos(a1-B1)-HD*sin(a1-B1)-W*cos(gamma);

F(13)=TD*sin(a1-B1)+HD*cos(a1-B1)-W*sin(gamma)+Rf*cos(gamma+tau)+Pc*sin(gamma+tau);
% F(10)=TD*sin(a1-B1)+HD*cos(a1-B1)-W*sin(gamma);

F(14)=Kh*(a1-B1)+W*h*sin(gamma)-Rf*h*cos(gamma+tau)+Mf-W*xcg*cos(gamma)-Rf*xcg*sin(gamma+tau)-...
    +Pc*l*cos(gamma+tau);
% F(11)=Kh*(a1-B1)+W*h*sin(gamma)-W*xcg*cos(gamma);
if strcmp(ptype,'autorot')
    F(15) = Q - rho*v_tip^2*A*R*sigma*(1/8*cd0*(1+mu^2)+...
    1/2*cla*(lambda*(1/3*theta0-1/2*lambda+mu*(1/2*a1-1/4*B1))));
end

end 