var x
    y
    z;

varexo u
       v;

parameters a1 a2 a3 a4
	   b1 b2 b3
	   c1;

a1 =  .50;
a2 =  .00;
a3 =  .70;
a4 =  .40;
b1 =  .90;
b2 =  .00;
b3 =  .80;
c1 =  .95;


model;
   y = a1*x(-1) + a2*x(1) + a3*z + a4*y(-1);
   z = b1*z(-1) + b2*z(1) + b3*x + u;
   x = c1*x(-1) + v;
end;

shocks;
var u; stderr .05;
var v; stderr .05;
end;

steady;

check;

stoch_simul(irf=0, periods=10000);

datatomfile('linear_model_data');


estimated_params;
a1,  .50;
a2,  .00;
a3,  .70;
a4,  .40;
b1,  .90;
b2,  .00;
b3,  .80;
c1,  .95;
stderr u, .05;
stderr v, .05;
end;

varobs y, z;


options_.particle.status = 1;
options_.particle.algorithm = 'sequential_importance_particle_filter';
options_.particle.initialization = 1;
options_.particle.pruning = 0;
options_.particle.number_of_particles = 20000;
options_.particle.resampling.status = 'systematic';
options_.particle.resampling.neff_threshold = .1;

options_.gstep(1) = 1e-4;
options_.gstep(2) = .1;

options_.mode_check_neighbourhood_size = 0.05;

set_dynare_threads('local_state_space_iteration_2',3);


estimation(order=2,nobs=100,datafile=linear_model_data,mh_replic=0);
