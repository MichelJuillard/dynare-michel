var Efficiency, efficiency;

varexo EfficiencyInnovation;

parameters rho, effstar, sigma2;

/*
** Calibration
*/


rho     =  0.950;
effstar =  1.000;
sigma2  =  0.0001;

external_function(name=mean_preserving_spread);

model;

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho));

end;

shocks;
var EfficiencyInnovation = sigma2;
end;

steady;

options_.ep.verbosity = 0;
options_.ep.stochastic.order = 0;
options_.ep.stochastic.nodes = 0;
options_.console_mode = 0;

ts = extended_path([],100);

options_.ep.verbosity = 0;
options_.ep.stochastic.order = 1;
options_.ep.stochastic.nodes = 3;
options_.console_mode = 0;

sts = extended_path([],100);

if max(max(abs(ts-sts)))>options_.dynatol.x
   disp('Stochastic Extended Path:: Something is wrong here (potential bug in extended_path.m)!!!')
end