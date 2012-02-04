var Capital, Output, Labour, Consumption, Efficiency, efficiency, ExpectedTerm;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, psi, delta, rho, effstar, sigma2;

/*
** Calibration
*/


beta    =  0.990;
theta   =  0.357;
tau     =  30.000;
alpha   =  0.450;
psi     =  -5.000;
delta   =  0.020;
rho     =  0.950;
effstar =  1.000;
sigma2  =  0.0001;

external_function(name=mean_preserving_spread);

model(use_dll);

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho));

  // Eq. n°3:
  Output = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour^psi))^(1/psi);

  // Eq. n°4:
  Consumption + Capital - Output - (1-delta)*Capital(-1);

  // Eq. n°5:
  ((1-theta)/theta)*(Consumption/(1-Labour)) - (1-alpha)*(Output/Labour)^(1-psi);

  // Eq. n°6:
  (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption - ExpectedTerm(1);

  // Eq. n°7:
  ExpectedTerm = beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital(-1))^(1-psi))+1-delta);

end;

shocks;
var EfficiencyInnovation = sigma2;
end;

steady;

options_.ep.verbosity = 0;
options_.console_mode = 0;

ts = extended_path([],1000);