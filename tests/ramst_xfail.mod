var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval(all_values_required);
c = aa*k^alph-delt*k;
end;

steady;

check;

shocks;
var x;
periods 1;
values 1.2;
end;

simul(periods=200);

rplot c;
rplot k;
