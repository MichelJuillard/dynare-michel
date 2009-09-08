var y, x;
varexo e;

parameters theta, rho, bet, xbar;

xbar = 0.0179;
rho = -0.139;
theta = -10;
bet = 0.95;

model;
y = bet*exp(theta*x(+1))*(1+y(+1));
x = (1-rho)*xbar + rho*x(-1) + e;
end;

initval;
x = 0.0179;
y = 0.3;
e = 0;
end;

vcov = [ 0.0012110];

order = 6;




