function [ys, info] = rbcii_steadystate(ys, exogenous)
% Steady state routine for rbc.mod (Business Cycle model with endogenous labour and CES production function)


% AUTHOR(S) 
%  stephane DOT adjemian AT univ DASH lemans DOT fr
%  frederic DOT karame AT univ DASH evry DOT fr    

% Output_per_unit_of_Capital = (((1/beta)-1+delta)/alpha)^(1/(1-psi));
% Consumption_per_unit_of_Capital = Output_per_unit_of_Capital - delta;
% Labour_per_unit_of_Capital = (((Output_per_unit_of_Capital/effstar)^psi-alpha)/(1-alpha))^(1/psi);
% Output_per_unit_of_Labour = Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
% Consumption_per_unit_of_Labour = Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;
% SteadyStateLabour = 1/(1 + Consumption_per_unit_of_Labour/((theta*(1-alpha)/(1-theta))*(Output_per_unit_of_Labour^(1-psi))));
% SteadyStateConsumption = Consumption_per_unit_of_Labour*SteadyStateLabour;
% SteadyStateCapital = SteadyStateLabour/Labour_per_unit_of_Capital;
% SteadyStateOutput =  Output_per_unit_of_Capital*SteadyStateCapital;
% ShareOfCapital = alpha/(alpha+(1-alpha)*Labour_per_unit_of_Capital^psi);
    
global M_
    
info = 0;

% Compute steady state ratios.
Output_per_unit_of_Capital=((1/M_.params(1)-1+M_.params(6))/M_.params(4))^(1/(1-M_.params(5)));
Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-M_.params(6);
Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/M_.params(8))^M_.params(5)-M_.params(4))/(1-M_.params(4)))^(1/M_.params(5));
Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

% Compute steady state share of capital.
ShareOfCapital=M_.params(4)/(M_.params(4)+(1-M_.params(4))*Labour_per_unit_of_Capital^M_.params(5));

% Compute steady state of the endogenous variables.
SteadyStateLabour=1/(1+Consumption_per_unit_of_Labour/((1-M_.params(4))*M_.params(2)/(1-M_.params(2))*Output_per_unit_of_Labour^(1-M_.params(5))));
SteadyStateConsumption=Consumption_per_unit_of_Labour*SteadyStateLabour;
SteadyStateCapital=SteadyStateLabour/Labour_per_unit_of_Capital;
SteadyStateOutput=Output_per_unit_of_Capital*SteadyStateCapital;

% Fill returned argument ys with steady state values.
ys = zeros(9,1);
ys(1)=SteadyStateCapital;
ys(2)=SteadyStateOutput;
ys(3)=SteadyStateLabour;
ys(4)=SteadyStateConsumption;
ys(5)=M_.params(8);
ys(7)=M_.params(1)*((((SteadyStateConsumption^M_.params(2))*((1-SteadyStateLabour)^(1-M_.params(2))))^(1-M_.params(3)))/SteadyStateConsumption)* ...
      (M_.params(4)*((SteadyStateOutput/SteadyStateCapital)^(1-M_.params(5)))+1-M_.params(6));