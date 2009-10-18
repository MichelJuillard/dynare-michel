// This file replicates the estimation of the CIA model from
// Frank Schorfheide (2000) "Loss function-based evaluation of DSGE models"
// Journal of  Applied Econometrics, 15, 645-670.
// the data are the ones provided on Schorfheide's web site with the programs.
// http://www.econ.upenn.edu/~schorf/programs/dsgesel.ZIP
// You need to have fsdat.m in the same directory as this file.
// This file replicates:
// -the posterior mode as computed by Frank's Gauss programs
// -the parameter mean posterior estimates reported in the paper
// -the model probability (harmonic mean) reported in the paper
// This file was tested with dyn_mat_test_0218.zip
// the smooth shocks are probably stil buggy
//
// The equations are taken from J. Nason and T. Cogley (1994)
// "Testing the implications of long-run neutrality for monetary business
// cycle models" Journal of Applied Econometrics, 9, S37-S70.
// Note that there is an initial minus sign missing in equation (A1), p. S63.
//
// Michel Juillard, February 2004
@#define bytecode = 1
@#define block = 1
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

@#if block == 1
  @#if bytecode == 1
    model(block, bytecode);
  @#else
    model(block);
  @#endif
@#else
  model;   
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
@#if bytecode == 1
  steady(solve_algo = 5);
@#else
  steady(solve_algo = 4);
@#endif

model_info;


@#if block == 0
  check;
@#endif

shocks;
var e_a;
periods 1;
values 0.16;
end;

@#if block == 1
  @#if bytecode == 1
    simul(periods=200, stack_solve_algo = 5);
  @#else
    simul(periods=200, stack_solve_algo = 1);
  @#endif
@#else
  simul(periods=200, stack_solve_algo = 0);
@#endif

rplot y;
rplot k;
rplot c;
