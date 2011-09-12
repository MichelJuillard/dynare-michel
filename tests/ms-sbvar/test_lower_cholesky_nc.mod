// same as test_upper_cholesky.mod, but with reordered variables. Results must be the same.
var R Pie Y;

varobs Y Pie R;

svar_identification;
lower_cholesky;
exclusion constants;
end;

sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4);

