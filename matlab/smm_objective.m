function [r,junk] = smm_objective(xparams,sample_moments,weighting_matrix,options,parallel)
% Evaluates the objective of the Simulated Moments Method.
%
% INPUTS:
%  xparams          [double]  p*1 vector of estimated parameters. 
%  sample_moments   [double]  n*1 vector of sample moments (n>=p).
%  weighting_matrix [double]  n*n symetric, positive definite matrix.
%  options          [      ]  Structure defining options for SMM.
%  parallel         [      ]  Structure defining the parallel mode settings (optional).
%
% OUTPUTS: 
%  r                [double]  scalar, the value of the objective function.
%  junk             [      ]  empty matrix.
%
% SPECIAL REQUIREMENTS
%  The user has to provide a file where the moment conditions are defined.

% Copyright (C) 2010 Dynare Team
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

global M_ options_
persistent mainStream mainState
persistent priorObjectiveValue
    
if nargin<5
    if isempty(mainStream)
        mainStream = RandStream.getDefaultStream;
        mainState  = mainStream.State;
    else
        mainStream.State = mainState;
    end
end

if isempty(priorObjectiveValue)
    priorObjectiveValue = Inf;
end

junk = [];

M_.params(options.estimated_parameters.idx) = xparams;

% Check for local determinacy of the deterministic steady state.
noprint = options_.noprint; options_.noprint = 1;
[local_determinacy_and_stability,info] = check; options_.noprint = noprint;
if ~local_determinacy_and_stability
    r = priorObjectiveValue * (1+info(2));
    return
end

simulated_moments = zeros(size(sample_moments));

if nargin<5
    for s = 1:options.number_of_simulated_sample
        time_series = extended_path([],options.simulated_sample_size,1);
        data = time_series(options.observed_variables_idx,options.burn_in_periods+1:options.simulated_sample_size);
        eval(['tmp = ' options.moments_file_name '(data);'])
        simulated_moments = simulated_moments + tmp;
        simulated_moments = simulated_moments / options.number_of_simulated_sample;
    end
else% parallel mode.
    job_number = 0;
    job_master = [];
    for i=1:length(parallel)
        job_remote = 0;
        for j=1:parallel(i).number_of_jobs
            job_remote = job_remote + 1; 
            job_number = job_number + 1;
            if strcmpi(hostname,remotename) && (job_number==1)
                master_job = [ job_number , i , job_remote ];
            else
                unix(['ssh -A ' parallel(i).login '@' parallel(i).machine './call_matlab_session.sh job' int2str(job_number) '&']);    
            end
        end
    end
    if ~isempty(master_job)
        tStartMasterJob = clock;
        eval(['job' int2str(master_job(1)) ';'])
        tElapsedMasterJob= etime(clock, tStartMasterJob);
    end
    tStart = clock;
    missing_jobs = 1;
    success_flag = 1;
    while missing_jobs
        tElapsed = etime(clock, tStart);
        if tElapsed>tElapsedMasterJob
            break
        end
        if ( length(dir('./intermediary_results_from_master_and_slaves/simulated_moments_slave_*.mat'))==job_number )
            missing_jobs = 0;
            success_flag = 0;
        end
    end
    if success_flag
        tmp = 0;
        for i=1:job_number
            load(['simulated_moments_slave_' int2str(i)]);
            tmp = tmp + simulated_moments;
        end
        simulated_moments = tmp/job_number;
    else
        r = priorObjectiveValue*1.1;
        return
    end
end

r = transpose(simulated_moments-sample_moments)*weighting_matrix*(simulated_moments-sample_moments);

priorObjectiveValue = r;