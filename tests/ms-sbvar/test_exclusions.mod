// same as test_lower_cholesky.mod, but using exclusion syntax
var R Pie Y;

varobs Y Pie R;

svar_identification;
exclusion lag 0;
equation 1, Pie, Y;
equation 2, Y;
end;

sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4);

