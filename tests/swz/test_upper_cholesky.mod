addpath '../../matlab/swz';
var Y Pie R;

model;
Y = 0;
Pie = 0;
R = 0;
end;

varobs Y Pie R;


sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4,restriction_fname=upper_cholesky);//for SBVAR

//ms_sbvar(datafile = data,freq=4,initial_year=1959,final_year=2005,nlags=4,restriction_fname=ftd_upperchol3v,
//         markov_file=specification_2v2c,mhm_file=MHM_input); //for markov switching