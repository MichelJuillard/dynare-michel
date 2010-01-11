function [param,sigma] = simulated_moments_estimation(xparam,dataset,options,parallel)
% Performs estimation by Simulated Moments Method.
%
% INPUTS:
%  xparam          [double]  p*1 vector of initial values for the estimated parameters. 
%  dataset         [      ]  Structure describing the data set.
%  options         [      ]  Structure defining options for SMM.
%  parallel        [      ]  Structure defining the parallel mode settings (optional). 
%
% OUTPUTS: 
%  param           [double]  p*1 vector of point estimates for the parameters.
%  sigma           [double]  p*p covariance matrix of the SMM estimates.
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

global M_ options_ oo_
    
% Load the dataset.
eval(dataset.name);
dataset.data = [];
for v = 1:dataset.number_of_observed_variables
    eval(['dataset.data = [ dataset.data , ' dataset.variables(v,:) ' ];'])
end
data = dataset.data(dataset.first_observation:dataset.first_observation+dataset.number_of_observations,:);

% Compute sample moments and the weighting matrix.
eval(['[sample_moments,weighting_matrix] = ' M_.fname '_moments;'])
weighting_matrix = inv(weighting_matrix);

% Initialize output.
sigma = [];
param = [];


% Set up parallel mode if needed.
if nargin>3
    if ~isunix
        error('The parallel version of SMM estimation is not implemented for non unix platforms!')
    end
    [junk,hostname] = unix('hostname --fqdn');    
    hostname = deblank(hostname);
    master_is_running_a_job = 0;
    for i=1:length(parallel)
        if strcmpi(hostname,parallel(i).machine)
            master_is_running_a_job = 1;
            break
        end
    end
    if ~master_is_running_a_job
        disp(' ')
        disp('Master has to run one job!')
        disp(' ')
    end
    disp(' ')
    disp('Master talks to its slaves...')
    disp(' ')
    % Save the workspace.
    save('master_variables.mat','options_','M_','oo_');
    % Send the workspace to each remote computer.
    disp('')
    for i = 1:length(parallel)
        if ~strcmpi(hostname,parallel(i).machine)
            unix(['scp master_variables.mat ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
        end
    end
    for i=1:length(parallel)
        % Write a bash script file to execute matlab scripts in the background
        fid = fopen('call_matlab_session.sh','w');
        fprintf(fid,'#!/bin/sh\n');
        fprintf(fid,'unset DISPLAY\n');
        fprintf(fid,['cd ' parallel(i).folder '\n']);
        fprintf(fid,['nohup ' parallel(i).matlab '/matlab -nodesktop -nodisplay -nojvm < $1 > /dev/null 2>&1\n']);
        fprintf(fid,'exit');
        fclose(fid);
        % Set the permission for this file (has to be executable)
        %fileattrib('call_matlab_session.sh','+x','u');
        unix(['chmod u+x call_matlab_session.sh']);
        % Send the script file on each remote computer
        if ~strcmpi(hostname,parallel(i).machine)
            unix(['scp call_matlab_session.sh ' , parallel(i).login , '@' , parallel(i).machine , ':~/' ]);
        else
            unix(['cp call_matlab_session.sh ~/call_matlab_session.sh']);
        end
    end
    % Send the files to each remote computer.
    for i = 1:length(parallel)
        if ~strcmpi(hostname,parallel(i).machine)
            unix(['scp ' M_.fname '_steadystate.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
            unix(['scp ' M_.fname '_moments.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
            unix(['scp ' M_.fname '_static.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
            if exist([M_.fname '_dynamic.c'])
                use_dll_flag = 1;
                unix(['scp ' M_.fname '_dynamic.c ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
            else
                use_dll_flag = 0;
                unix(['scp ' M_.fname '_dynamic.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
            end
        end
    end
    % If needed, compile dynamic model mex file on each remote computer
    if use_dll_flag
        % Write a matlab script that will trigger the compilation of the mex file.
        fid = fopen('compile_model.m', 'w');
        fprintf(fid,[' eval(''mex -O LDFLAGS=''''-pthread -shared -Wl,--no-undefined'''' ' M_.fname '_dynamic.c'')  ']);
        fprintf(fid, '\n exit');
        fclose(fid);
        for i = 1:length(parallel)
            if ~strcmpi(hostname,parallel(i).machine)
                % Send the generated matlab script to the remote computer.
                unix(['scp compile_model.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder]);
                % Compile the mex file on the remote computer.
                unix(['ssh ' parallel(i).login , '@' , parallel(i).machine ' ./call_matlab_session.sh  compile_model.m']);
            end
        end
    end
    % Write the matlab script files for the evaluation of the simulated moment conditions
    job = 0;
    for i=1:length(parallel)
        for j=1:parallel(i).number_of_jobs
            job = job+1;
            % Create random number streams
            write_job(hostname, parallel(i).machine, parallel(i).dynare, ...
                      options.simulated_sample_size, length(sample_moments), ...
                      dataset.observed_variables_idx, options.estimated_parameters.idx, options.burn_in_periods, [M_.fname '_moments'], parallel(i).number_of_simulations, ...
                      parallel(i).number_of_threads_per_job, job, j);
            if ~strcmpi(hostname,parallel(i).machine)
                unix(['scp ' , 'job' , int2str(job) , '.m ' , parallel(i).login , '@' , parallel(i).machine , ':' parallel(i).folder ]);
            end
        end
    end
    disp(' ')
    disp('... And slaves do as ordered.')
    disp(' ')
    if exist('intermediary_results_from_master_and_slaves','dir')
        unix('rm -rf intermediary_results_from_master_and_slaves');
    end
    unix('mkdir intermediary_results_from_master_and_slaves');
    unix('chmod -R u+x intermediary_results_from_master_and_slaves');
end

disp('');

if options.optimization_routine==1
    % Set options for csminwel. 
    H0 = 1e-4*eye(options.estimated_parameters.nb);
    ct = 1e-4;
    it = 1000;
    vb = 2;
    % Minimization of the objective function.
    if nargin==3
        [fval,param,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
            csminwel('smm_objective',xparam,H0,[],ct,it,2,options_.gradient_epsilon,sample_moments,weighting_matrix,options);    
    elseif nargin>3
        [fval,param,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
            csminwel('smm_objective',xparam,H0,[],ct,it,2,options_.gradient_epsilon,sample_moments,weighting_matrix,options,parallel);
    end
elseif options.optimization_routine==2
    optim_options = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-4,'TolX',1e-4);
    if isfield(options_,'optim_opt')
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if nargin==3
        [param,fval,exitflag] = fminsearch('smm_objective',xparam,optim_options,sample_moments,weighting_matrix,options);
    else
        [param,fval,exitflag] = fminsearch('smm_objective',xparam,optim_options,sample_moments,weighting_matrix,options,parallel);
    end
end


function write_job(hostname, remotename, dynare_path, sample_size, number_of_moments, observed_variables_idx, parameters_idx, burn_in_periods, moments_file_name, number_of_simulations,threads_per_job, slave_number, job_number)
 
fid = fopen(['job' int2str(slave_number) '.m'],'w');

fprintf(fid,['% Generated by ' hostname '.\n\n']);

if ( strcmpi(hostname,remotename) && (job_number>1) )  || ~strcmpi(hostname,remotename) 
    fprintf(fid,'load(''master_variables'');\n');
    fprintf(fid,'assignin(''base'',''M_'',M_);\n');
    fprintf(fid,'assignin(''base'',''oo_'',oo_);\n');
    fprintf(fid,'assignin(''base'',''options_'',options_);\n\n');
end

if ( strcmpi(hostname,remotename) && (job_number>1) ) || ~strcmpi(hostname,remotename)
    fprintf(fid,['addpath ' dynare_path '\n']);
    fprintf(fid,['dynare_config;\n\n']);
end

fprintf(fid,['simulated_moments = zeros(' int2str(number_of_moments) ',1);\n\n']);

fprintf(fid,'load(''estimated_parameters.mat'');\n');
fprintf(fid,['M_.params([' num2str(parameters_idx)  ']) = xparams;\n\n']);

fprintf(fid,['stream=RandStream(''mt19937ar'',''Seed'',' int2str(slave_number) ');\n']);
fprintf(fid,['RandStream.setDefaultStream(stream);\n\n']);

fprintf(fid,['maxNumCompThreads(' int2str(threads_per_job) ');\n\n']);

fprintf(fid,['for s = 1:' int2str(number_of_simulations) '\n'] );
fprintf(fid,['    time_series = extended_path([],' int2str(sample_size) ',1);\n']);
fprintf(fid,['    data = time_series([' int2str(observed_variables_idx) '],' int2str(burn_in_periods) '+1:' int2str(sample_size) ');\n']);
fprintf(fid,['    eval(''tmp = ' moments_file_name '(data);'');\n']);
fprintf(fid,['    simulated_moments = simulated_moments + tmp;\n']);
fprintf(fid,['end;\n\n']);

fprintf(fid,['simulated_moments = simulated_moments/' int2str(number_of_simulations) ';\n']);
fprintf(fid,['save(''simulated_moments_slave_' int2str(slave_number) '.dat'',''simulated_moments'',''-ascii'');\n']);

if ~strcmpi(hostname,remotename)
    fprintf(fid,['unix(''scp simulated_moments_slave_' int2str(slave_number) '.dat ' hostname ':' pwd '/intermediary_results_from_master_and_slaves '');\n']);
    fprintf(fid,['unix(''rm simulated_moments_slave_' int2str(slave_number) '.dat'');\n']);
else
    fprintf(fid,['unix(''cp simulated_moments_slave_' int2str(slave_number) '.dat '  'intermediary_results_from_master_and_slaves '');\n']);
    fprintf(fid,['unix(''rm simulated_moments_slave_' int2str(slave_number) '.dat'');\n']);
end

if ((job_number>1) && strcmpi(hostname,remotename)) || ~strcmpi(hostname,remotename) 
    fprintf(fid,'exit');
end

fclose(fid);