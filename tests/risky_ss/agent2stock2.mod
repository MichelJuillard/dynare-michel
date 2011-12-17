var c1 c2 x1 x2 p1 p2 d1 d2 y1 y2;
varexo eps_1 eps_2 eta_1 eta_2;

parameters beta, gamma, kappa, rho_y, rho_d;

beta = 0.96;
gamma = 4;
kappa = -0.5;
rho_y = 0.9;
rho_d = 0.9;

model;
p1*c1^(-gamma-1) = beta*c1(+1)^(-gamma-1)*(d1(+1)+p1(+1));
p2*c1^(-gamma-1) = beta*c1(+1)^(-gamma-1)*(d2(+1)+p2(+1));
p1*c2^(-gamma-1) = beta*c2(+1)^(-gamma-1)*(d1(+1)+p1(+1));
p2*c2^(-gamma-1) = beta*c2(+1)^(-gamma-1)*(d2(+1)+p2(+1));
c1 = y1 - (x1(-1)*(p1+d1)-x1*p1) + (x2(-1)*(p2+d2)-x2*p2); 		
c2 = y2 + (x1(-1)*(p1+d1)-x1*p1) - (x2(-1)*(p2+d2)-x2*p2);
y1 = (1-rho_y)*0.5 + rho_y*y1(-1) + eps_1;
y2 = (1-rho_y)*0.5 + rho_y*y2(-1) + eps_2;
d1 = (1-rho_d)*0.5 + rho_d*d1(-1) + eta_1;
d2 = (1-rho_d)*0.5 + rho_d*d2(-1) + eta_2;
end;

shocks;
var eps_1; stderr 0.01;
var eps_2; stderr 0.01;
var eta_1; stderr 0.01;
var eta_2; stderr 0.01;
corr eps_1,eta_1 = -0.5;
corr eps_2,eta_2 = -0.5;
end;

initval;
c1 = 0.5;
c2 = 0.5;
y1 = 0.5;
y2 = 0.5;
d1 = 0.5;
d2 = 0.5;
p1 = 12.51;
p2 = 12.51;
x1 = 0.26;
x2 = 0.25;
end;

options_.risky_steadystate = 1;

stoch_simul(irf=0);
