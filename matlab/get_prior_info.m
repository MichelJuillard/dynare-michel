function get_prior_info(info)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
    global options_ M_ estim_params_ oo_
    
    if ~nargin
        info = 0;
    end
    
    [xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);
    plot_priors(bayestopt_,M_,options_);
    
    PriorNames = { 'Beta' , 'Gamma' , 'Gaussian' , 'Inverted Gamma' , 'Uniform' , 'Inverted Gamma -- 2' };
    
    if size(M_.param_names,1)==size(M_.param_names_tex,1)% All the parameters have a TeX name.
        fidTeX = fopen('priors_data.tex','w+');
        % Column 1: a string for the name of the prior distribution. 
        % Column 2: the prior mean.
        % Column 3: the prior standard deviation.
        % Column 4: the lower bound of the prior density support.
        % Column 5: the upper bound of the prior density support.
        % Column 6: the lower bound of the interval containing 80% of the prior mass. 
        % Column 7: the upper bound of the interval containing 80% of the prior mass.
        prior_trunc_backup = options_.prior_trunc ;
        options_.prior_trunc = (1-options_.prior_interval)/2 ;
        PriorIntervals = prior_bounds(bayestopt_) ;
        options_.prior_trunc = prior_trunc_backup ;
        for i=1:size(bayestopt_.name,1)
            [tmp,TexName] = get_the_name(i,1);
            PriorShape = PriorNames{ bayestopt_.pshape(i) };
            PriorMean = bayestopt_.pmean(i);
            PriorStandardDeviation = bayestopt_.pstdev(i);
            switch bayestopt_.pshape(i)
              case { 1 , 5 }
                LowerBound = bayestopt_.p3(i);
                UpperBound = bayestopt_.p4(i);
              case { 2 , 4 , 6 }
                LowerBound = bayestopt_.p3(i);
                UpperBound = '$\infty$';
              case 3
                if isinf(bayestopt_.p3(i))
                    LowerBound = '$-\infty$';
                else
                    LowerBound = bayestopt_.p3(i);
                end
                if isinf(bayestopt_.p4(i))
                    UpperBound = '$\infty$';
                else
                    UpperBound = bayestopt_.p4(i);
                end
              otherwise
                error('get_prior_info:: Dynare bug!')
            end
            format_string = build_format_string(bayestopt_,i);
            fprintf(fidTeX,format_string, ...
                    TexName, ...
                    PriorShape, ...
                    PriorMean, ...
                    PriorStandardDeviation, ...
                    LowerBound, ...
                    UpperBound, ...
                    PriorIntervals(i,1), ...
                    PriorIntervals(i,2) );
        end
        fclose(fidTeX);
    end
    
    M_.dname = M_.fname;
    
    if info% Prior simulations.
       results = prior_sampler(1,M_,bayestopt_,options_,oo_);
       results.prior.mass
    end
    
    
    
    
    
    
function format_string = build_format_string(bayestopt,i)
    format_string = ['%s & %s & %6.4f &'];
    if isinf(bayestopt.pstdev(i))
        format_string = [ format_string , ' %s &'];
    else
        format_string = [ format_string , ' %6.4f &'];
    end
    if isinf(bayestopt.p3(i))
        format_string = [ format_string , ' %s &'];
    else
        format_string = [ format_string , ' %6.4f &'];
    end
    if isinf(bayestopt.p4(i))
        format_string = [ format_string , ' %s &'];
    else
        format_string = [ format_string , ' %6.4f &'];
    end
    format_string = [ format_string , ' %6.4f & %6.4f \\\\ \n'];