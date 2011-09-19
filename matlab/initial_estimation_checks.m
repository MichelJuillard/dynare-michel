function DynareResults = initial_estimation_checks(xparam1,DynareDataset,Model,EstimatedParameters,DynareOptions,BayesInfo,DynareResults)
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

if DynareDataset.info.nvobs>Model.exo_nbr+EstimatedParameters.nvn
    error(['initial_estimation_checks:: Estimation can''t take place because there are less shocks than observed variables!'])
end

if DynareOptions.dsge_var
    [fval,cost_flag,info] = DsgeVarLikelihood(xparam1,gend);
else
    [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults);
end

% when their is an analytical steadystate, check that the values
% returned by *_steadystate match with the static model
if DynareOptions.steadystate_flag
    [ys,check] = feval([Model.fname '_steadystate'],...
                       DynareResults.steady_state,...
                       [DynareResults.exo_steady_state; ...
                        DynareResults.exo_det_steady_state]);
    if size(ys,1) < Model.endo_nbr
        if length(Model.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,Model.aux_vars,...
                                                        Model.fname,...
                                                        DynareResults.exo_steady_state,...
                                                        DynareResults.exo_det_steady_state,...
                                                        Model.params,...
                                                        DynareOptions.bytecode);
        else
            error([Model.fname '_steadystate.m doesn''t match the model']);
        end
    end
    DynareResults.steady_state = ys;
    % Check if the steady state obtained from the _steadystate file is a
    % steady state.
    check1 = 0;
    if isfield(DynareOptions,'unit_root_vars') && DynareOptions.diffuse_filter == 0
        if isempty(DynareOptions.unit_root_vars)
            if ~DynareOptions.bytecode
                check1 = max(abs(feval([Model.fname '_static'],...
                                       DynareResults.steady_state,...
                                       [DynareResults.exo_steady_state; ...
                                    DynareResults.exo_det_steady_state], Model.params))) > DynareOptions.dynatol ;
            else
                [info, res] = bytecode('static','evaluate',DynareResults.steady_state,...
                                       [DynareResults.exo_steady_state; ...
                                    DynareResults.exo_det_steady_state], Model.params);
                check1 = max(abs(res)) > DynareOptions.dynatol;
            end
            if check1
                error(['The seadystate values returned by ' Model.fname ...
                       '_steadystate.m don''t solve the static model!' ])
            end
        end
    end
end

if info(1) > 0
    disp('Error in computing likelihood for initial parameter values')
    print_info(info, DynareOptions.noprint)
end

if any(abs(DynareResults.steady_state(BayesInfo.mfys))>1e-9) && (DynareOptions.prefilter==1)
    disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
    disp(['variables using demeaned data!'])
    error('You should change something in your mod file...')
end

disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
