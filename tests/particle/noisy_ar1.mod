var x y;
varexo u v;

parameters c, r;

c = .500;
r = .975;

model;
y = x + u;
x = r*x(-1) + 0.0000001*x(1) + v;
end;

shocks;
var u; stderr sqrt(0.02);
var v; stderr sqrt(0.02);
end;

steady;

check;

stoch_simul(irf=0, periods=10000);

datatomfile('noisy_ar1_data');

estimated_params;
c,  .500;
r,  .975;
stderr u, sqrt(.2);
stderr v, sqrt(.02);
end;

varobs y;


particle = 1;

if particle;
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

else;

estimation(nobs=100,datafile=noisy_ar1_data,mh_replic=0);

end;
