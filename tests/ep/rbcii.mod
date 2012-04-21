@#define extended_path_version = 1

var Capital, Output, Labour, Consumption, Investment, Efficiency, efficiency, ExpectedTerm;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, psi, delta, rho, effstar, sigma2;

/*
** Calibration
*/


beta    =  0.990;
theta   =  0.357;
tau     =  2.000;
alpha   =  0.450;
psi     =  -0.500;
delta   =  0.020;
rho     =  0.995;
effstar =  1.000;
sigma2  =  0.001;


@#if extended_path_version
    rho = 0.800;
@#endif

external_function(name=mean_preserving_spread);

model(use_dll);

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho));

  // Eq. n°3:
  Output = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour^psi))^(1/psi);

  // Eq. n°4:
  Capital = max(Output-Consumption,0) + (1-delta)*Capital(-1);

  // Eq. n°5:
  ((1-theta)/theta)*(Consumption/(1-Labour)) - (1-alpha)*(Output/Labour)^(1-psi);

  // Eq. n°6:
  ExpectedTerm = beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital(-1))^(1-psi))+1-delta);

  // Eq. n°7:
  Investment = Capital - (1-delta)*Capital(-1);

  // Eq. n°8: (Euler equation, to be skipped if investment is on its lower bound)
  (Investment>0)*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption - ExpectedTerm(1)) + (1-(Investment>0))*(Output-Consumption);

end;


@#if extended_path_version

    shocks;
    var EfficiencyInnovation = sigma2;
    end;

    steady;

    options_.maxit_ = 100;
    options_.ep.verbosity = 0;
    options_.ep.stochastic.order = 0;
    options_.ep.stochastic.nodes = 2;
    options_.console_mode = 0;
    ts = extended_path([],100);

    options_.ep.stochastic.order = 1;
    sts = extended_path([],100);

    figure(1)
    plot(ts(2,:)-ts(4,:));

    figure(2)
    plot(sts(2,:)-sts(4,:));

    figure(3)
    plot(sts(2,:)-ts(2,:))

@#else

    shocks;
    var EfficiencyInnovation;
    periods 1;
    values -.4;
    end;

    steady;

    options_.maxit_ = 100;

    simul(periods=4000);

    n = 100;

    figure('Name','(rbcii) Investment.');
    plot(Output(1:n)-Consumption(1:n),'-b','linewidth',2)

    figure('Name','(rbcii) Lagrange multiplier associated to the positivity constraint on investment.');
    plot(LagrangeMultiplier(1:n),'-b','linewidth',2)

@#endif