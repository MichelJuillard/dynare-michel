function [i_var,nvar] = varlist_indices(varlist)
% function [i_var,nvar] = varlist_indices(varlist)
% returns the indices of a list of endogenous variables
%
% INPUT
%   varlist:    (character area) list of variables
%
% OUTPUT
%   i_var:      variable indices in M_.endo_names
%   nvar:       number of variables in varlist
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2004-2008)
% Gnu Public License.

    global M_
    
    endo_nbr = M_.endo_nbr;
    
    if isempty(varlist)
        i_var = (1:endo_nbr)';
        nvar = endo_nbr;
    else
        i_var = [];
        for i=1:size(varlist,1)
            tmp = strmatch(varlist(i,:),M_.endo_names,'exact');
            if isempty(tmp)
                error([tmp ' isn''t an endogenous variable'])
            end
            i_var = [i_var; tmp];
        end
        nvar = length(i_var);
    end
    