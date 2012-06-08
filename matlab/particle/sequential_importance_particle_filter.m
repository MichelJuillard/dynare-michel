function [LIK,lik] = sequential_importance_particle_filter(ReducedForm,Y,start,DynareOptions)
% Evaluates the likelihood of a nonlinear model with a particle filter (optionally with resampling).

%@info:
%! @deftypefn {Function File} {@var{y}, @var{y_} =} sequential_importance_particle_filter (@var{ReducedForm},@var{Y}, @var{start}, @var{DynareOptions})
%! @anchor{particle/sequential_importance_particle_filter}
%! @sp 1
%! Evaluates the likelihood of a nonlinear model with a particle filter (optionally with resampling).
%!
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ReducedForm
%! Structure describing the state space model (built in @ref{non_linear_dsge_likelihood}).
%! @item Y
%! p*smpl matrix of doubles (p is the number of observed variables), the (detrended) data.
%! @item start
%! Integer scalar, likelihood evaluation starts at observation 'start'.
%! @item DynareOptions
%! Structure specifying Dynare's options.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item LIK
%! double scalar, value of (minus) the logged likelihood.
%! @item lik
%! smpl*1 vector of doubles, density of the observations at each period.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @ref{non_linear_dsge_likelihood}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2012 Dynare Team
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

% AUTHOR(S) frederic DOT karame AT univ DASH evry DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

persistent init_flag
persistent mf0 mf1
persistent number_of_particles
persistent sample_size number_of_observed_variables number_of_structural_innovations

% Set default value for start
if isempty(start)
    start = 1;
end

% Get steady state and mean.
steadystate = ReducedForm.steadystate;
constant = ReducedForm.constant;
state_variables_steady_state = ReducedForm.state_variables_steady_state;

% Set persistent variables (if needed).
if isempty(init_flag)
    mf0 = ReducedForm.mf0;
    mf1 = ReducedForm.mf1;
    sample_size = size(Y,2);
    number_of_observed_variables = length(mf1);
    number_of_structural_innovations = length(ReducedForm.Q);
    number_of_particles = DynareOptions.particle.number_of_particles;
    init_flag = 1;
end

% Set local state space model (first order approximation).
ghx  = ReducedForm.ghx;
ghu  = ReducedForm.ghu;

% Set local state space model (second order approximation).
ghxx = ReducedForm.ghxx;
ghuu = ReducedForm.ghuu;
ghxu = ReducedForm.ghxu;

% Get covariance matrices.
Q = ReducedForm.Q;
H = ReducedForm.H;
if isempty(H)
    H = 0;
end

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);
Q_lower_triangular_cholesky = chol(Q)';

% Set seed for randn().
set_dynare_seed('default');

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables;
lik  = NaN(sample_size,1);

% Initialization of the weights across particles.
nb_obs_resamp = 0 ;
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean);
for t=1:sample_size
    yhat = bsxfun(@minus,StateVectors,state_variables_steady_state);
    epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);
    PredictedObservedMean = mean(tmp(mf1,:),2);
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    dPredictedObservedMean = bsxfun(@minus,tmp(mf1,:),PredictedObservedMean);
    PredictedObservedVariance = (dPredictedObservedMean*dPredictedObservedMean')/number_of_particles+H;
    lnw = -.5*(const_lik+log(det(PredictedObservedVariance))+sum(PredictionError.*(PredictedObservedVariance\PredictionError),1));
    dfac = max(lnw);
    wtilde = weights.*exp(lnw-dfac);
    lik(t) = log(sum(wtilde))+dfac;
    weights = wtilde/sum(wtilde);
    if strcmp(DynareOptions.particle.resampling.status,'generic')
        Neff = 1/(weights*weights');
    end
    if (strcmp(DynareOptions.particle.resampling.status,'generic') && Neff<DynareOptions.particle.resampling.neff_threshold*sample_size ) || strcmp(DynareOptions.particle.resampling.status,'systematic')
        nb_obs_resamp = nb_obs_resamp+1 ;
        StateVectors = tmp(mf0,resample(weights,DynareOptions.particle.resampling.method1,DynareOptions.particle.resampling.method2));
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif strcmp(DynareOptions.particle.resampling.status,'smoothed')
        StateVectors = multivariate_smooth_resampling(weights',tmp(mf0,:)',number_of_particles,DynareOptions.particle.resampling.number_of_partitions)';
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif strcmp(DynareOptions.particle.resampling.status,'none')
        StateVectors = tmp(mf0,:);
    end
end

LIK = -sum(lik(start:end));