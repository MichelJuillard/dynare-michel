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

model(block,bytecode);

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
options_.ep.stochastic = 0;
options_.console_mode = 0;


ts = extended_path([],1000);


options_.ep.verbosity = 0;
options_.ep.stochastic = 1;
options_.console_mode = 0;

sts = extended_path([],1000);


if max(max(abs(ts-sts)))>1e-12
   disp('Stochastic Extended Path:: Something is wrong here (potential bug in extended_path.m)!!!')
end