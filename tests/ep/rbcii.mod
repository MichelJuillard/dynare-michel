@#define extended_path_version = 1

var Capital, Output, Labour, Consumption,  Investment, Output1, Labour1, Consumption1, Output2, Labour2, Consumption2, Efficiency, efficiency, ExpectedTerm;

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

external_function(name=mean_preserving_spread,nargs=2);

model;

  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho,sigma2));

  (((Consumption1^theta)*((1-Labour1)^(1-theta)))^(1-tau))/Consumption1 - ExpectedTerm(1);

  ExpectedTerm = beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital(-1))^(1-psi))+1-delta);

  ((1-theta)/theta)*(Consumption1/(1-Labour1)) - (1-alpha)*(Output1/Labour1)^(1-psi);

  Output1 = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour1^psi))^(1/psi);

  Consumption2  = Output2;

  ((1-theta)/theta)*(Consumption2/(1-Labour2)) - (1-alpha)*(Output2/Labour2)^(1-psi);

  Output2 = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour2^psi))^(1/psi);

  Consumption = (Output1 > Consumption1)*Consumption1 + (1-(Output1 > Consumption1))*Consumption2;

  Labour = (Output1 > Consumption1)*Labour1 + (1-(Output1 > Consumption1))*Labour2;

  Output = (Output1 > Consumption1)*Output1 + (1-(Output1 > Consumption1))*Output2;

  Capital = Output-Consumption + (1-delta)*Capital(-1);

  Investment = Capital - (1-delta)*Capital(-1);

end;

// Write analytical steady state file (without globals)
options_.steadystate_flag = 2;
copyfile('rbcii_steady_state.m','rbcii_steadystate2.m');

@#if extended_path_version

    shocks;
    var EfficiencyInnovation = 1;
    end;

    steady(nocheck);

    options_.maxit_ = 100;
    options_.ep.verbosity = 0;
    options_.ep.stochastic.order = 0;
    options_.ep.stochastic.nodes = 2;
    options_.console_mode = 0;
    ts = extended_path([],100);

    options_.ep.stochastic.order = 1;
    sts = extended_path([],100);

//    options_.ep.stochastic.order = 2;
//    sts2 = extended_path([],100);

//    options_.ep.stochastic.order = 3;
//    sts3 = extended_path([],100);

//    figure(1)
//    plot(ts(2,:)-ts(4,:));

//    figure(2)
//    plot(sts(2,:)-sts(4,:));

//    figure(3)
//    plot(sts(2,:)-ts(2,:))

//    figure(4)
//    plot([(ts(2,:)-ts(4,:))' (sts(2,:)-sts(4,:))' (sts2(2,:)-sts2(4,:))' (sts3(2,:)-sts3(4,:))']) 
//    plot([(ts(2,:)-ts(4,:))' (sts(2,:)-sts(4,:))']) 

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