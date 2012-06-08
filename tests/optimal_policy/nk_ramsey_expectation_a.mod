//MODEL:
// test for using expectation operator with Ramsey policy
// Here, the expected variable is replaced manually with an auxiliary variable (epai)
// We check that the result is the same as with an explicit expectation operator in nk_ramsey_expectation.mod
// Note that the example doesn't make any sense from an economic point of view.

//------------------------------------------------------------------------------------------------------------------------
//1. Variable declaration
//------------------------------------------------------------------------------------------------------------------------

var pai, c, n, r, a, epai;

//4 variables + 1 shock

varexo u;




//------------------------------------------------------------------------------------------------------------------------
// 2. Parameter declaration and calibration
//-------------------------------------------------------------------------------------------------------------------------

parameters beta, rho, epsilon, omega, phi, gamma;

beta=0.99;
gamma=3;       //Frish elasticity
omega=17;        //price stickyness
epsilon=8;      //elasticity for each variety of consumption
phi=1;           //coefficient associated to labor effort disutility

rho=0.95;  //coefficient associated to productivity shock


//-----------------------------------------------------------------------------------------------------------------------
// 3. The model
//-----------------------------------------------------------------------------------------------------------------------


model;


a=rho*(a(-1))+u;

1/c=beta*(1/(c(+1)))*(r/epai);               //euler


omega*pai*(pai-1)=beta*omega*(c/(c(+1)))*(pai(+1))*(pai(+1)-1)+epsilon*exp(a)*n*(c/exp(a)*phi*n^gamma-(epsilon-1)/epsilon);  //NK pc
//pai*(pai-1)/c = beta*pai(+1)*(pai(+1)-1)/c(+1)+epsilon*phi*n^(gamma+1)/omega-exp(a)*n*(epsilon-1)/(omega*c);  //NK pc

(exp(a))*n=c+(omega/2)*((pai-1)^2);

epai = pai(+1);
end;

//--------------------------------------------------------------------------------------------------------------------------
// 4. Steady state
//---------------------------------------------------------------------------------------------------------------------------

initval;
% this is the exact steady state under optimal policy
%important for the comparison

pai=1;
epai=1;
r=1/beta;
c=((epsilon-1)/(epsilon*phi))^(1/(1+gamma));
n=c;
a=0;


end;



//---------------------------------------------------------------------------------------------------------------------------
// 5. shocks
//---------------------------------------------------------------------------------------------------------------------------

shocks;
var u; stderr 0.008;

end;

//--------------------------------------------------------------------------------------------------------------------------
// 6. Ramsey problem
//--------------------------------------------------------------------------------------------------------------------------

planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));

options_.solve_tolf=1e-12;
ramsey_policy(planner_discount=0.99);

o1=load('nk_ramsey_expectation_results');
if (norm(o1.oo_.dr.ghx-oo_.dr.ghx,inf) > 1e-12)
   error('ghx doesn''t match')
end
if (norm(o1.oo_.dr.ghu-oo_.dr.ghu,inf) > 1e-12)
   error('ghu doesn''t match')
end
if (abs(o1.oo_.planner_objective_value(1)-oo_.planner_objective_value(1)) > 1e-12)
   error('planner objective value doesn''t match')
end

