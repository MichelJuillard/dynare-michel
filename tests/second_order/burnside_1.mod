var y x;
varexo e;

parameters beta theta rho xbar;
xbar = 0.0179;
rho = -0.139;
theta = -1.5;
theta = -10;
beta = 0.95;

model;
y = beta*exp(theta*x(+1))*(1+y(+1));
x = (1-rho)*xbar + rho*x(-1)+e;
end;

shocks;
var e; stderr 0.0348;
end;

initval;
x = xbar;
y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

resid(1);

steady;

stoch_simul(order=2,irf=0);

sigma2=M_.Sigma_e;
i = linspace(1,800,800);
a = theta*xbar*i+(theta^2*sigma2)/(2*(1-rho)^2)*(i-2*rho*(1-rho.^i)/(1-rho)+rho^2*(1-rho.^(2*i))/(1-rho^2));
a1 = theta^2*sigma2/(2*(1-rho)^2)*(i-2*rho*(1-rho.^i)/(1-rho)+rho^2*(1-rho.^(2*i))/(1-rho^2));
b = theta*rho*(1-rho.^i)/(1-rho);

dr = oo_.dr;
x1 = [dr.ghx(2); dr.ghu(2); dr.ghxx(2); dr.ghuu(2); dr.ghxu(2); dr.ghs2(2)];
x2 = [
 sum(beta.^i.*exp(theta*xbar*i).*b*rho)
 sum(beta.^i.*exp(theta*xbar*i).*b)
 sum(beta.^i.*exp(theta*xbar*i).*(b*rho).^2)
 sum(beta.^i.*exp(theta*xbar*i).*b.^2)
 sum(beta.^i.*exp(theta*xbar*i).*b.^2*rho)
 sum(beta.^i.*exp(theta*xbar*i).*a1*2)
];
if any(abs(x1-x2) > 1e-14);
   error('burnside_1 doesn''t reproduce the analytical solution');
end;