function print_info(info)

% function print_info(info)
% Prints error messages
%
% INPUTS
%    info=1:         the model doesn't determine the current variables '...' uniquely
%    info=2:         MJDGGES returns the following error code'
%    info=3:         Blanchard Kahn conditions are not satisfied: no stable '...' equilibrium
%    info=4:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy
%    info=5:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy due to rank failure
%    info=20:        Impossible to find the steady state. Either the model' ...' doesn't have 
%                    a unique steady state of the guess values' ...' are too far from the solution
%    info=30:        Variance can't be computed
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2005-2007)
% Gnu Public License.


  global options_

  options_ = set_default_option(options_,'noprint',0);
  
  if ~options_.noprint
    if info(1) > 10 & info(1) < 20
      disp('Failure in dr_algo=2')
      info(1) = info(1) - 10;
    end
    switch info(1)
     case 1
      error(['The model doesn''t determine the current variables' ...
	     ' uniquely'])
     case 2
      error(['MJDGGES returns the following error code' ...
	     int2str(info(2))])
     case 3
      error(['Blanchard Kahn conditions are not satisfied: no stable' ...
	     ' equilibrium'])
     case 4
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy'])
     case 5
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy due to rank failure'])
     case 20
      error(['Impossible to find the steady state. Either the model' ...
	     ' doesn''t have a unique steady state of the guess values' ...
	     ' are too far from the solution']) 
     case 30
      error('Variance can''t be computed')
     otherwise
      error('This case shouldn''t happen. Contact the authors of Dynare')
    end
  end
  