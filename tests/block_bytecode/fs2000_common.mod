var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

@#if block_bytecode == 2
model(block, bytecode);
@#else
@# if block_bytecode == 1
model(block);
@# else
model;
@# endif
@#endif

/*0*/  exp(gam+e_a) = dA ;
/*1*/  log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
/*2*/  -P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
/*3*/  l/n = W;
/*4*/  -(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
/*5*/  R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
/*6*/  1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
/*7*/  c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
/*8*/  P*c = m; 
/*9*/  m-1+d = l;
/*10*/ e = exp(e_a);
/*11*/ k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a)) = y ;
/*12*/ gy_obs = dA*y/y(-1);
/*13*/ gp_obs = (P/P(-1))*m(-1)/dA;
end;

initval;
k = 6;
m = mst;
P = 2.25;
c = 0.45;
e = 1;
W = 4;
R = 1.02;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gy_obs = exp(gam);
gp_obs = exp(-gam);
dA = exp(gam);
e_a=0;
e_m=0;
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

options_.maxit_=10;
steady(solve_algo = @{solve_algo});

@#if block_bytecode > 0
model_info;
@#endif

shocks;
var e_a;
periods 1;
values 0.16;
end;

simul(periods=200, stack_solve_algo = @{stack_solve_algo});

@#if block_bytecode > 0
if ~exist('fs2000_simk_results.mat','file');
   error('fs2000_simk must be run first');
end;

oo1 = load('fs2000_simk_results','oo_');

err = max(max(abs(oo_.endo_simul - oo1.oo_.endo_simul)))
disp(['Max error in simulation: ' num2str(err)])
if err > 1e-4
   error('Error above the threshold')
end;
@#endif
