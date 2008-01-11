function [xparams, logpost] = GetOneDraw(type)

% function [xparams, logpost] = GetOneDraw(type)
% draws one row from metropolis
%
% INPUTS
%    type:      posterior
%               prior
%        
% OUTPUTS
%    xparams:   vector of estimated parameters (drawn from posterior distribution)
%    logpost:   log of the posterior density relative to this row
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2005-2008)
% Gnu Public License.

  switch type
   case 'posterior'
    [xparams, logpost] = metropolis_draw(0);
   case 'prior'
    xparams = prior_draw(0);
  end