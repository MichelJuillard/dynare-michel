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

shocks;
var e_a; stderr 0.035;
end;

steady;

//stoch_simul(order=2,drop=0,periods=250,noprint,nograph) y l i ;
//disp([y l i ]) ;
//disp(oo_.mean) ;

estimated_params;
alp, uniform_pdf,,, 0.0001, 0.5;
bet, uniform_pdf,,, 0.0001, 0.99;
tet, uniform_pdf,,, 0.0001, 1;
tau, uniform_pdf,,, 0.0001, 100;
delt, uniform_pdf,,, 0.0001, 0.05;
rho, uniform_pdf,,, 0.8, 0.99;
stderr e_a, uniform_pdf,,, 0.00001, 0.1;
stderr y, uniform_pdf,,, 0.00001, 0.1;
stderr l, uniform_pdf,,, 0.00001, 0.1;
stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

estimated_params_init;
alp, 0.4;
bet, 0.97;
tet, 0.357 ;
tau, 50;
delt, 0.02;
rho, 0.9 ;
stderr e_a, .035;
stderr y, .0175;//.00158;
stderr l, .00312;//.0011;
stderr i, .00465;//.000866;
end;

varobs y l i ;

//options_.gstep(1) = 1e-4;
//options_.gstep(2) = .1;

options_.particle.status = 1;
options_.particle.initialization = 1;
options_.particle.pruning = 0;
options_.particle.number_of_particles = 1000 ;
options_.particle.resampling.status = 'systematic';

//options_.particle.resampling.method1 = 'traditional' ;
//options_.particle.resampling.method1 = 'residual' ;
options_.particle.resampling.method1 = 'smooth' ;

options_.particle.reampling.method2 = 'kitagawa' ;
//options_.particle.resampling.method2 = 'stratified' ;

options_.particle.resampling.number_of_partitions = 1;
options_.particle.resampling.neff_threshold = .5;

set_dynare_threads('local_state_space_iteration_2',3);

options_.particle.algorithm = 'sequential_importance_particle_filter';
//options_.particle.algorithm = 'auxiliary_particle_filter';
//options_.particle.algorithm = 'gaussian_mixture_filter';
//options_.particle.algorithm = 'conditional_particle_filter';
//options_.particle.algorithm = 'gaussian_filter';

//options_.particle.IS_approximation_method = 'quadrature' ;
options_.particle.IS_approximation_method = 'cubature' ;
//options_.particle.IS_approximation_method = 'unscented' ;

//options_.particle.approximation_method = 'quadrature' ;
options_.particle.approximation_method = 'cubature' ;
//options_.particle.approximation_method = 'unscented' ;
//options_.particle.approximation_method = 'MonteCarlo' ;

//options_.mh_posterior_mode_estimation=1 ;

// online
options_.particle.liu_west_delta = 0.9 ;
options_.mode_check_node_number = 250 ;

estimation(datafile=data_risky_perturb3,nograph,order=2,nobs=100,mh_replic=0,mode_compute=7,mode_check);
