% Copyright (C) 2001 Michel Juillard
function result = check
global M_ oo_ options_ it_

tempex = oo_.exo_simul;
if ~options_.initval_file
	oo_.exo_simul = ones(M_.maximum_lead+M_.maximum_lag+1,1)*transpose(oo_.exo_steady_state);
end

% check if ys is steady state
it_ = M_.maximum_lag+1;
fh = str2func([M_.fname '_static']);
if max(abs(feval(fh,oo_.steady_state,oo_.exo_simul))) > options_.dynatol
	if exist([M_.fname '_steadystate'])
    	[dr.ys, cheik] = feval([M_.fname '_steadystate'],oo_.steady_state,oo_.exo_simul);
	else
    	[dr.ys, cheik] = dynare_solve([M_.fname '_static'], ...
		    oo_.steady_state,oo_.exo_simul);
	end
	if cheik
    	error('CHECK: convergence problem in DYNARE_SOLVE')
	end
else
	dr.ys = oo_.steady_state;
end

dr = dr1(1,dr,1);
oo_.exo_simul = tempex;
    
nyf = nnz(dr.kstate(:,2)>M_.maximum_lag+1);
[m_lambda,i]=sort(abs(oo_.eigenvalues));
disp(' ')
disp('EIGENVALUES:')
disp(sprintf('%16s %16s %16s\n','Modulus','Real','Imaginary'))
z=[m_lambda real(oo_.eigenvalues(i)) imag(oo_.eigenvalues(i))]';
disp(sprintf('%16.4g %16.4g %16.4g\n',z))
options_ = set_default_option(options_,'qz_criterium',1.000001);
disp(sprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', ...
      nnz(abs(oo_.eigenvalues) > options_.qz_criterium)));
disp(sprintf('for %d forward-looking variable(s)',nyf));
disp(' ')
if dr.rank == nyf
	disp('The rank condition is verified.')
else
	disp('The rank conditions ISN''T verified!')
end
disp(' ')
  
  % keep oo_.eigenvalues for backward compatibility
  %  oo_.eigenvalues = oo_.eigenvalues;

  % 2/9/99 MJ: line 15, added test for absence of exogenous variable.
  % 8/27/2000 MJ: change JACOB call. Added ...,1 to cumsum()
  % 6/24/01 MJ: added count of abs(eigenvalues) > 1
  % 2/21/02 MJ: count eigenvalues > 1[+1e-5]
  % 01/22/03 MJ: warning(warning_state) needs parentheses for Matlab 6.5
  % 03/20/03 MJ: changed name of global from lambda to oo_.eigenvalues to avoid
  %              name conflicts with parameters
  % 05/21/03 MJ: replace computation by dr1.m and add rank check
  % 06/05/03 MJ: corrected bug when M_.maximum_lag > 0
  % 22/01/05 SA: variable check --> cheik.
