function [ivar,vartan] = set_stationary_variables_list()
% This function builds a vector of indices targeting to the stationary
% variables in varlist.
% 
% INPUTS 
%   None.
%  
% OUTPUTS 
%   o ivar       [integer]  nvar*1 vector of indices (nvar is the number
%                           of stationary variables).
%   o vartan     [char]     array of characters (with nvar rows).
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright Dynare Team (2007)
% Gnu Public License.    
global options_ M_    
varlist = options_.varlist;
if isempty(varlist)
    varlist = options_.varobs;
    options_.varlist = varlist;
end
nvar = rows(varlist);
if ~isempty(options_.unit_root_vars)
    vartan = [];
    for i=1:nvar
        if isempty(strmatch(deblank(varlist(i,:)),options_.unit_root_vars,'exact'))       
            vartan = strvcat(vartan,varlist(i,:));
        end
    end
else
    vartan = varlist;
end
nvar = size(vartan,1);
ivar = zeros(nvar,1);
for i = 1:nvar
    ivar(i) = strmatch(deblank(vartan(i,:)),M_.endo_names,'exact');
end