var R Pie Y;

model;
Y = 0;
Pie = 0;
R = 0;
end;

varobs Y Pie R;

svar_identification;
lower_cholesky;
end;

markov_switching(chain=1,number_of_states=2,duration=2.5);

svar(variances, chain=1);

ms_estimation(datafile=data,freq=4,initial_year=1959,final_year=2005,nlags=4);
ms_simulation(mh_replic=1000);
ms_compute_mdd;
ms_compute_probabilities;
