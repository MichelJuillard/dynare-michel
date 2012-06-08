// Test preprocessor interface for extended path

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

external_function(name=mean_preserving_spread,nargs=2);

model(use_dll);

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho,sigma2));

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

steady_state_model;
efficiency = 0;
Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho,sigma2));
// Compute steady state ratios.
Output_per_unit_of_Capital=((1/beta-1+delta)/alpha)^(1/(1-psi));
Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-delta;
Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/Efficiency)^psi-alpha)/(1-alpha))^(1/psi);
Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

// Compute steady state share of capital.
ShareOfCapital=alpha/(alpha+(1-alpha)*Labour_per_unit_of_Capital^psi);

/// Compute steady state of the endogenous variables.
Labour=1/(1+Consumption_per_unit_of_Labour/((1-alpha)*theta/(1-theta)*Output_per_unit_of_Labour^(1-psi)));
Consumption = Consumption_per_unit_of_Labour*Labour;
Capital = Labour/Labour_per_unit_of_Capital;
Output = Output_per_unit_of_Capital*Capital;
ExpectedTerm = beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital)^(1-psi))+1-delta);
end;


shocks;
var EfficiencyInnovation = sigma2;
end;

steady;

extended_path(periods=100);

