// This file deals with the resolution and estimation of a basic DSGE model with
//employment for comparison with the benchmark in Gauss which solves with
//the same particular filter but global methodology.
//
// January 2010

var k A c l i y;
varexo e_a;

parameters alp bet tet tau delt rho ;
alp = 0.4;
bet = 0.99;
tet = 0.357 ;
tau = 50 ;
delt = 0.02;
rho = 0.95;

model;
c = ((1 - alp)*tet/(1-tet))*A*(1-l)*((k(-1)/l)^alp) ;
y = A*(k(-1)^alp)*(l^(1-alp)) ;
i = y-c ;
k = (1-delt)*k(-1) + i ;
log(A) = rho*log(A(-1)) + e_a ;
(((c^(tet))*((1-l)^(1-tet)))^(1-tau))/c - bet*((((c(+1)^(tet))*((1-l(+1))^(1-tet)))^(1-tau))/c(+1))*(1 -delt+alp*(A(1)*(k^alp)*(l(1)^(1-alp)))/k)=0 ;
end;

steady;
resid(1);


shocks;
var e_a; stderr 0.007;
end;

steady;


estimated_params;
alp, uniform_pdf,,, 0.0001, 1;
bet, uniform_pdf,,, 0.75, 0.999;
tet, uniform_pdf,,, 0.0001, 1;
tau, uniform_pdf,,, 0.0001, 100;
delt, uniform_pdf,,, 0.0001, 0.05;
rho, uniform_pdf,,, 0.0001, 0.999;
stderr e_a, uniform_pdf,,, 0.00001, 0.1;
stderr y, uniform_pdf,,, 0.00001, 0.1;
stderr l, uniform_pdf,,, 0.00001, 0.1;
stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

estimated_params_init;
alp, 0.4;
bet, 0.99;
tet, 0.357 ;
tau, 50;
delt, 0.02;
rho, 0.95;
stderr e_a, .035;
stderr y, .0175;//.00158;
stderr l, .00312;//.0011;
stderr i, .00465;//.000866;
end;


varobs y l i;

options_.particle.status = 1;
options_.particle.algorithm = 'sequential_importance_particle_filter';
options_.particle.initialization = 1;
options_.particle.pruning = 1;
options_.particle.number_of_particles = 20000;
options_.particle.resampling.status = 'systematic';
options_.particle.resampling.neff_threshold = .1;

options_.gstep(1) = 1e-4;
options_.gstep(2) = .1;

options_.mode_check_neighbourhood_size = 0.05;

set_dynare_threads('local_state_space_iteration_2',3);

estimation(datafile=data_risky_perturb2,nograph,order=2,nobs=100,mh_replic=0,mode_compute=7,mode_check);
