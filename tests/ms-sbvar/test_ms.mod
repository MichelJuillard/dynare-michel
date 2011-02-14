// same as test_upper_cholesky.mod, but with reordered variables. Results must be the same.
var R Pie Y;

model;
Y = 0;
Pie = 0;
R = 0;
end;

varobs Y Pie R;


sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4,restriction_fname=lower_cholesky);

