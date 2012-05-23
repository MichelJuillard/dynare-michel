var c k aa;
varexo x;

parameters alph gam delt bet rho;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
rho=0.9;

model;
c + k - aa*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa(+1)*alph*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
log(aa) = rho*log(aa(-1))+x;
end;

initval;
x = 0;
aa = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady(transform_logpow);

check;

shocks;
var x;
periods 1;
values 0.2;
end;

simul(periods=200);

rplot c;
rplot k;

write_latex_dynamic_model;
write_latex_static_model;