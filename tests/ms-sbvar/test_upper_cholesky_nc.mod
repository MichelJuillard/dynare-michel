var Y Pie R;

model;
Y = 0;
Pie = 0;
R = 0;
end;

varobs Y Pie R;

svar_identification;
exclusion constants;
upper_cholesky;
end;

sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4);

