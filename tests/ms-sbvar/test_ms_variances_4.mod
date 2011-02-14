// same as test_upper_cholesky.mod, but with reordered variables. Results must be the same.
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

markov_switching(chain=1,number_of_states=4,duration=2.5);

svar(variances, chain=1);

ms_sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4,draws_nbr_modified_harmonic_mean=10000);

