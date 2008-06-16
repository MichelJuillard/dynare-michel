function moments=compute_model_moments(dr,options)

% function compute_model_moments(options)
% Computes posterior filter smoother and forecasts
%
% INPUTS
%    dr:        structure describing model solution
%    options:   structure of Dynare options
%    
% OUTPUTS
%    moments: a cell array containing requested moments
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

% subset of variables
    varlist = options.varlist;

    moments = th_autocovariances(dr,varlist);
