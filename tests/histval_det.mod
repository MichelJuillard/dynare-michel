// Tests the histval block in deterministic setup
// In particular test if it works on endos and exos substituted by an aux var

var c k z;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
z = 0.9*z(-2);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(-1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
z = 0;
end;

steady;

check;

histval;
k(0) = 0.6;
x(0) = 0.9;
z(-1) = 0.1;
end;

simul(periods=200);
