function sol = extractdata(ptype,x,input)
sol.TD=x(1); sol.HD=x(2); sol.Rf=x(3); sol.Mf=x(4); sol.Pc=x(5);
sol.a1=x(6); sol.gamma=x(7); sol.mu=x(8); sol.lambda=x(9);
sol.alfaD=x(10); sol.u=x(11); sol.v1=x(12);

if strcmp(ptype,'inverse')
    sol.theta0=x(13); sol.B1=x(14);
    sol.V = input.V ; sol.tau = input.tau;
elseif strcmp(ptype,'direct')
    sol.V=x(13); sol.tau=x(14);
    sol.B1 = input.B1; sol.theta0 = input.theta0;
    a1_H = sol.a1-sol.B1; 
    sol.TH = sol.TD*cos(a1_H)-sol.HD*sin(a1_H);
    sol.HH = sol.TD*sin(a1_H)+sol.HD*cos(a1_H); 
    sol.a1_H = a1_H*180/pi;
else
    sol.theta0=x(13); sol.B1=x(14); sol.tau = x(15);
    sol.V = input.V ; sol.Q = input.Q;
end

end