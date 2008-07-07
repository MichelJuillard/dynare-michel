function options=set_default_option(options,field,default)

% function options=set_default_option(options,field,default)
% Sets the option value 
% 
% INPUTS
%    options
%    field:   option name
%    default: assigns a value
%    
% OUTPUTS
%    options
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  if ~isfield(options,field)
    options.(field) = default;
  end
  
  % 06/07/03 MJ added ; to eval expression