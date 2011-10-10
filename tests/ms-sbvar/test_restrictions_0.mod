// same as test_lower_cholesky.mod, but using exclusion syntax
var R Pie Y;

varobs Y Pie R;

svar_identification;
restriction equation 1, coeff(Pie,0) = 0;
restriction equation 1, coeff(Y,0) = 0;
restriction equation 2, coeff(Y,0) = 0;
end;

sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4);

