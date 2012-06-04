var y;
varexo eps;
parameters rho1 rho2;

rho1 = 0.1;

rho2 = 0.2;

model;
log(y) = rho1*log(y(-1)) + rho2*log(y(-2)) + eps;
end;

initval;
y=1;
eps=0;
end;

steady;
check;

shocks;
var eps;
periods 1:9;
values -0.0104;
end;

simul(periods=100);
