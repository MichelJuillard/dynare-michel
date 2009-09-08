// this is a file trying to replicate example1.mod as optimization of a social planner
// it serves also as an example combining +2 lead and optimal policy

var Y, C, K, A, H, B;

varexo EPS, NU;

parameters rho, beta, alpha, delta, theta, psi, tau;
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 1/(1.03^0.25);
delta = 0.025;
psi   = 0;
theta = 2.95;

planner_objective log(C)-theta*H^(1+psi)/(1+psi);

planner_discount beta;

model;
//Y = exp(A)*K^alpha*H^(1-alpha);
Y = exp(A)*exp(A(+1))*exp(A(+2))*K^alpha*H^(1-alpha);
K = exp(B(-1))*(Y(-1)-C(-1)) + (1-delta)*K(-1);
A = rho*A(-1) + tau*B(-1) + EPS;
B = tau*A(-1) + rho*B(-1) + NU;
end;

initval;
A = 0;
B = 0;
H = ((1-alpha)/(theta*(1-(delta*alpha)/(1/beta-1+delta))))^(1/(1+psi));
Y = (alpha/(1/beta-1+delta))^(alpha/(1-alpha))*H;
K = alpha/(1/beta-1+delta)*Y;
C = Y - delta*K;
end;

vcov = [
  0.0002  0.00005;
  0.00005 0.0001
];

order = 2;
