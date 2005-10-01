function initial_estimation_checks(xparam1,gend,data)
global dr1_test bayestopt_ estim_params_ options_ oo_ M_

nv = size(data,1);
if nv-size(options_.varobs,1)
	disp(' ')
	disp(['Declared number of observed variables = ' int2str(size(options_.varobs,1))])
	disp(['Number of variables in the database   = ' int2str(nv)])
    disp(' ')
	error(['Estimation can''t take place because the declared number of observed' ...
	  'variables doesn''t match the number of variables in the database.'])
end
if nv > M_.exo_nbr+estim_params_.nvn
    error(['Estimation can''t take place because there are less shocks than' ...
	  'observed variables'])
end
r = rank(data);
if r < nv
    error(['Estimation can''t take place because the data are perfectly' ...
 	   ' correlated']);
end

fval = DsgeLikelihood(xparam1,gend,data);
if exist(dr1_test)
	disp(dr1_test)
	switch(dr1_test(1))
		case 1
			error('The steady state can''t be found');
		case 2
			error(['Estimation can''t take place because there are an infinity of' ...
					' stable solutions']);
		case 3
			error(['Estimation can''t take place because there is no stable' ...
					' solution']);
		case 4
			error(['Estimation can''t take place because of singularity in Kalman' ...
					' filter']);
		otherwise
	end
end
