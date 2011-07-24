function adaptive_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin)
%function adaptive_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin)
% Random walk Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o varargin              list of argument following mh_bounds
%  
% OUTPUTS 
%   None  
%
% ALGORITHM 
%   Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code sutable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle),
% is removed from this function and placed in random_walk_metropolis_hastings_core.m funtion.
% Then the DYNARE parallel package contain a set of pairs matlab functions that can be executed in
% parallel and called name_function.m and name_function_core.m. 
% In addition in parallel package we have second set of functions used
% to manage the parallel computation.
%
% This function was the first function to be parallelized, later other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management funtions.

% Copyright (C) 2006-2011 Dynare Team
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

global M_ options_ bayestopt_ estim_params_ oo_


old_options = options_;

accept_target = options_.amh.accept_target;

options_.mh_jscale = tune_scale_parameter(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin{:});

for i=1:options_.amh.cova_steps
    options_.mh_replic = options_.amh.cova_replic;
    random_walk_metropolis_hastings(TargetFun,ProposalFun, ...
                                    xparam1,vv,mh_bounds,varargin{:});
    CutSample(M_,options_,estim_params_);
    [junk,vv] = compute_mh_covariance_matrix();
    bayestopt_.jscale = tune_scale_parameter(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin{:});
end

options_.mh_replic = old_options.mh_replic;
record = random_walk_metropolis_hastings(TargetFun,ProposalFun, ...
                                         xparam1,vv,mh_bounds,varargin{:});
                


function selected_scale = tune_scale_parameter(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin)
global options_ bayestopt_

selected_scale = [];

maxit = options_.amh.scale_tuning_maxit;
accept_target = options_.amh.accept_target;
test_runs = options_.amh.scale_tuning_test_runs;
tolerance = options_.amh.scale_tuning_tolerance;
Scales = zeros(maxit,1);
AvRates = zeros(maxit,1);
Scales(1) = options_.mh_jscale;

for i=1:maxit
    options_.mh_replic = options_.amh.scale_tuning_blocksize;
    bayestopt_.jscale = Scales(i);
    record = random_walk_metropolis_hastings(TargetFun,ProposalFun, ...
                                             xparam1,vv, ...
                                             mh_bounds,varargin{:});
    AvRates(i) = mean(record.AcceptationRates);    
    disp(AvRates(1:i)')

    if i >= test_runs
        a_mean = mean(AvRates((i-test_runs+1):i));
        if (a_mean > (1-tolerance)*accept_target) && ...
                (a_mean < (1+tolerance)*accept_target)
            selected_scale = mean(Scales((i-test_runs+1):i));
            disp(['Selected scale: ' num2str(selected_scale)])
            return
        end
    end
    
    if i == 1
        if AvRates(1) > accept_target
            Scales(i+1) = 2*Scales(i);
        else
            Scales(i+1) = Scales(i)/2;
        end
    elseif i < maxit
        X = [ones(i,1) Scales(1:i)];
        b = X\(AvRates(1:i)-accept_target);
        Scales(i+1) = -b(1)/b(2);
        if Scales(i+1) < 0.05
            Scales(i+1) = 0.05;
        elseif Scales(i+1) > 2
            Scales(i+1) = 2;
        end
    else
        error('AMH scale tuning: tuning didn''t converge')
    end
    
    options_.load_mh_file = 1;
    
    disp(Scales(1:i)')
end