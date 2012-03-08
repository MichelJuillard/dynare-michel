var x y z;

varexo u v;

parameters a d e f g;

a = 0.98;
g = 0.25;
d = 0.80;
e = 0.90;
f = 0.50;

model;

    z = a*z(-1) + u + g*u(-1);

    y = d*y(1) + z;

    x = e*x(-1) + f*y + v;

end;

steady;

check;

shocks;
var u = 0.01;
var v = 0.01;
end;

stoch_simul(order=1,irf=0,periods=10000);

save('mydata.mat','x','y','z');

estimated_params;
  e, beta_pdf, 0.90, 0.05;
  d, beta_pdf, 0.80, 0.05;
  g, beta_pdf, 0.25, 0.05;
end;

varobs x, y;

//estimation(datafile=mydata,order=1,first_obs=5001,nobs=100,mh_replic=0);

options_.particle.status = 1;
options_.particle.algorithm = 'sequential_importance_particle_filter';
options_.particle.initialization = 1;
options_.particle.number_of_particles = 1000;

options_.mode_check_neighbourhood_size = 0.05;

set_dynare_threads('local_state_space_iteration2',2);

estimated_params_init;
e, 0.9009;
d, 0.7912;
g, 0.2448;
end;

estimation(datafile=mydata,order=2,first_obs=5001,nobs=100,mh_replic=0,mode_compute=8);
