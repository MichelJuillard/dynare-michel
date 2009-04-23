var c k a;
varexo e;

parameters alpha delta beta rho;

alpha = 0.33;
delta = 0.025;
beta = 0.99;
rho = 0.9;

model;
1/c = beta*(1/c(+1))*(a(+1)*alpha*k^(alpha-1)+1-delta);
c+k = a*k(-1)^alpha + (1-delta)*k(-1);
log(a) = rho*log(a(-1))+e;
end;

initval;
c = 2;
k = 28;
a = 1;
end;


steady;

shocks;
var e; stderr 0.01;
end;

options_.use_k_order=1;
stoch_simul(irf=0);
