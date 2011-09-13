var R Pie Y;

varobs Y Pie R;

svar_identification;
lower_cholesky;
end;

markov_switching(chain=1,number_of_states=2,duration=2.5);

svar(variances, chain=1);

set_dynare_seed(5);

ms_estimation(datafile=data
		,freq=4
		,initial_year=1959
		,final_year=2005
		,nlags=4
		,max_repeated_optimization_runs=1
		,max_number_of_stages=0
                ,lower_cholesky
);
ms_simulation(mh_replic=1000);
ms_compute_mdd;
ms_compute_probabilities;
ms_irf;
ms_forecast;
ms_variance_decomposition;
