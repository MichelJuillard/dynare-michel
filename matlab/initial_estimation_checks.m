function oo = initial_estimation_checks(xparam1,dataset,M,estim_params,options,bayestopt,oo); %initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations
% function initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
% Checks data (complex values, ML evaluation, initial values, BK conditions,..)
%
% INPUTS
%    xparam1: vector of parameters to be estimated
%    gend:    scalar specifying the number of observations
%    data:    matrix of data
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2011 Dynare Team
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

if dataset.info.nvobs>M.exo_nbr+estim_params.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because there are less shocks than observed variables!'])
end

if options.dsge_var
    [fval,cost_flag,info] = DsgeVarLikelihood(xparam1,gend);
else
    [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
end

% when their is an analytical steadystate, check that the values
% returned by *_steadystate match with the static model
if options.steadystate_flag
    [ys,check] = feval([M.fname '_steadystate'],...
                       oo.steady_state,...
                       [oo.exo_steady_state; ...
                        oo.exo_det_steady_state]);
    if size(ys,1) < M.endo_nbr
        if length(M.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M.aux_vars,...
                                                        M.fname,...
                                                        oo.exo_steady_state,...
                                                        oo.exo_det_steady_state,...
                                                        M.params,...
                                                        options.bytecode);
        else
            error([M.fname '_steadystate.m doesn''t match the model']);
        end
    end
    oo.steady_state = ys;
    % Check if the steady state obtained from the _steadystate file is a
    % steady state.
    check1 = 0;
    if isfield(options,'unit_root_vars') && options.diffuse_filter == 0
        if isempty(options.unit_root_vars)
            if ~options.bytecode
                check1 = max(abs(feval([M.fname '_static'],...
                                       oo.steady_state,...
                                       [oo.exo_steady_state; ...
                                    oo.exo_det_steady_state], M.params))) > options.dynatol ;
            else
                [info, res] = bytecode('static','evaluate',oo.steady_state,...
                                       [oo.exo_steady_state; ...
                                    oo.exo_det_steady_state], M.params);
                check1 = max(abs(res)) > options.dynatol;
            end
            if check1
                error(['The seadystate values returned by ' M.fname ...
                       '_steadystate.m don''t solve the static model!' ])
            end
        end
    end
end

if info(1) > 0
    disp('Error in computing likelihood for initial parameter values')
    print_info(info, options.noprint)
end

if any(abs(oo.steady_state(bayestopt.mfys))>1e-9) && (options.prefilter==1)
    disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
    disp(['variables using demeaned data!'])
    error('You should change something in your mod file...')
end

disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
