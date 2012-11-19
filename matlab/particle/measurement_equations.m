function measure = measurement_equations(StateVectors,ReducedForm,DynareOptions) 
mf1 = ReducedForm.mf1;
ghx  = ReducedForm.ghx(mf1,:);
ghu  = ReducedForm.ghu(mf1,:);
ghxx = ReducedForm.ghxx(mf1,:);
ghuu = ReducedForm.ghuu(mf1,:);
ghxu = ReducedForm.ghxu(mf1,:);
constant = ReducedForm.constant(mf1,:);
state_variables_steady_state = ReducedForm.state_variables_steady_state;
number_of_structural_innovations = length(ReducedForm.Q);
yhat = bsxfun(@minus,StateVectors,state_variables_steady_state) ;
measure = local_state_space_iteration_2(yhat,zeros(number_of_structural_innovations,size(yhat,2)),ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);



