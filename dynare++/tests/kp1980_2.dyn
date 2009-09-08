// Model from Kydland & Prescott JEDC 1980

// case 2: time inconsistent optimal policy with different policy and consumer objectives

var C G K TAU Z;

varexo EPS;

parameters eta beta alpha delta phi a rho; 

eta = 2;
beta = 0.99;
alpha = 0.3;
delta = 0.10;
phi = 2.5;
a = 0.1;
rho = 0.7;

planner_objective C^(1-eta)/(1-eta) + a*G^(1-phi)/(1-phi);

planner_discount beta;

model;
K = (1-delta)*K(-1) + (exp(Z(-1))*K(-1)^alpha - C(-1) - G(-1));
G = TAU*alpha*K^alpha;
Z = rho*Z(-1) + EPS;
C^(-eta) = beta*C(+1)^(-eta)*(1-delta+exp(Z(+1))*alpha*K(+1)^(alpha-1)*(1-alpha*TAU(+1)));
end;

initval;
TAU = 0.70;
K = ((delta+1/beta-1)/(alpha*(1-alpha*TAU)))^(1/(alpha-1));
G = TAU*alpha*K^alpha;
C =  K^alpha - delta*K - G;
Z = 0;
end;

order = 4;

vcov = [
	0.01
];
