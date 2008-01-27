function x = get_param_by_name(pname)

% function x = get_param_by_name(pname)
% returns the value of a parameter identified by its name
%  
% INPUTS:
%   pname:  parameter name
%
% OUTPUTS
%   x:      parameter value
%
% SPECIAL REQUIREMENTS
%   none
%
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License  

  global M_
  
  i = strmatch(pname,M_.param_names,'exact');
  
  if isempty(i)
    error(['Can''t find parameter ' pname])
  end
  
  x = M_.params(i);