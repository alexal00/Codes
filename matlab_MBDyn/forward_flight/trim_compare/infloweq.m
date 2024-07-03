function F = infloweq(x,data,cT)
lambda = x(1);
mu = data.mu*cos(data.alfaD);
alfaD = data.alfaD;
F(1) = - lambda + mu*tan(alfaD)+ cT/(2*sqrt(mu^2+lambda^2));

end