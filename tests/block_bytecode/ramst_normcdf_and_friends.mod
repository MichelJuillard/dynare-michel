// Tests the normcdf(), normpdf() and erf() functions, in the static and the dynamic byecode file

var c k t u v;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model(bytecode, block);
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
t = normcdf(x, 2, 3);
u = normpdf(x, 1, 0.5);
v = erf(x);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
t = 0;
u = 0;
v = 0;
end;

steady(solve_algo=5);

//check;

shocks;
var x;
periods 1;
values 1.2;
end;

simul(periods=20, stack_solve_algo=5);

if(abs(oo_.steady_state(5) - erf(1)) > 1e-10)
   error('Test failed in static bytecode for erf')
end

if (abs(oo_.endo_simul(5, 2) - erf(1.2)) > 1e-10)
   error('Test failed in dynamic bytecode for erf')
end

if(abs(oo_.steady_state(4) - normpdf(1, 1, 0.5)) > 1e-10)
   error('Test failed in static bytecode for normpdf')
end

if (abs(oo_.endo_simul(4, 2) - normpdf(1.2, 1, 0.5)) > 1e-10)
   error('Test failed in dynamic bytecode for normpdf')
end

if (abs(oo_.steady_state(3) - normcdf(1, 2, 3)) > 1e-10)
   error('Test failed in static bytecode for normcdf')
end

if (abs(oo_.endo_simul(3, 2) - normcdf(1.2, 2, 3)) > 1e-10)
   error('Test failed in dynamic bytecode for normcdf')
end

