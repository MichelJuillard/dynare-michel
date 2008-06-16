function s=get_moments_size(options)

% function PosteriorFilterSmootherAndForecast(Y,gend, type)
% Computes posterior filter smoother and forecasts
%
% INPUTS
%    options: structure of Dynare options
%    
% OUTPUTS
%    s: size of moments for a given model and options
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

    global M_
    
    n = size(options.varlist,1);
    
    if n == 0
        n = M_.endo_nbr;
    end
    
    n2 = n*n;

    s = n; % mean
    s = s + n;  % std errors
    s = s + n2; % variance
    s = s + n2; % correlations
    s = s + options.ar*n2; % auto-correlations