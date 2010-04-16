// Tests the normcdf() function, in the static M-file, and in a dynamic C file

var c k t u v w;
varexo x;

parameters alph gam delt bet aa c_steady_state;
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
w = steady_state(k);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
t = 0;
u = 0;
v = 0;
w = 0;
end;

steady(solve_algo=5);

//check;

shocks;
var x;
periods 1;
values 1.2;
end;

simul(periods=20, stack_solve_algo=5);

if(abs(oo_.steady_state(2) - oo_.endo_simul(6,2)) > 1e-10)
   error('Test failed in bytecode for steady_state')
end

