// Example 1 from Collard's guide to Dynare
var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_1 = oo_.endo_simul;

set_dynare_seed('reset')
stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_2 = oo_.endo_simul;

set_dynare_seed('default')
stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_3 = oo_.endo_simul;

t1 = endo_simul_1-endo_simul_2;
t2 = endo_simul_1-endo_simul_3;

if any(abs(t1(:))>1e-12) || any(abs(t2(:))>1e-12)
   disp('Test 1.')
   error('Test failure:: Problem with the seed of the random number algorithm')
end

set_dynare_seed(57)
stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_4 = oo_.endo_simul;

set_dynare_seed('reset')
stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_5 = oo_.endo_simul;

t3 = endo_simul_4-endo_simul_5;

if any(abs(t3(:))>1e-12)
   disp('Test 2.')
   error('Test failure:: Problem with the seed of the random number algorithm')
end

set_dynare_seed('default')
stoch_simul(periods=1000,irf=0,nomoments);
endo_simul_6 = oo_.endo_simul;

t4 = endo_simul_6-endo_simul_1;

if any(abs(t4(:))>1e-12)
   disp('Test 3.')
   error('Test failure:: Problem with the seed of the random number algorithm')
end

if ~exist('OCTAVE_VERSION') && ~matlab_ver_less_than(7.7)

    set_dynare_seed('mlfg6331_64',0)
    stoch_simul(periods=1000,irf=0,nomoments);
    endo_simul_7 = oo_.endo_simul;

    set_dynare_seed('reset')
    stoch_simul(periods=1000,irf=0,nomoments);
    endo_simul_8 = oo_.endo_simul;

    t5 = endo_simul_7-endo_simul_8;

    if any(abs(t5(:))>1e-12)
        disp('Test 4.')
        error('Test failure:: Problem with the seed of the random number algorithm')
    end

    set_dynare_seed('default')
    stoch_simul(periods=1000,irf=0,nomoments);
    endo_simul_9 = oo_.endo_simul;
    
    t6 = endo_simul_9-endo_simul_1;

    if any(abs(t6(:))>1e-12)
        disp('Test 5.')
        error('Test failure:: Problem with the seed of the random number algorithm')
    end

end
