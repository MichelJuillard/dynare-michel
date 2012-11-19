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

% Copyright (C) 2011, 2012 Dynare Team
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

% AUTHOR(S) frederic DOT karame AT univ DASH lemans DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

persistent init_flag
persistent mf0 mf1
persistent number_of_particles
persistent sample_size number_of_observed_variables number_of_structural_innovations

% Set default value for start
if isempty(start)
    start = 1;
end

% Set flag for prunning
pruning = DynareOptions.particle.pruning;

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

% Initialization of the likelihood.
const_lik = log(2*pi)*number_of_observed_variables;
lik  = NaN(sample_size,1);

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
if pruning
    StateVectorMean_ = StateVectorMean;
    StateVectorVarianceSquareRoot_ = StateVectorVarianceSquareRoot;
end

% Get the rank of StateVectorVarianceSquareRoot
state_variance_rank = size(StateVectorVarianceSquareRoot,2);

% Factorize the covariance matrix of the structural innovations
Q_lower_triangular_cholesky = chol(Q)';
[PredictedStateMean,PredictedStateVarianceSquareRoot,StateVectorMean,StateVectorVarianceSquareRoot] = ...
        gaussian_filter_bank(ReducedForm,Y(:,1),StateVectorMean,StateVectorVarianceSquareRoot,Q_lower_triangular_cholesky,Q_lower_triangular_cholesky,H,DynareOptions) ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean) ;

% Set seed for randn().
set_dynare_seed('default');

% Initialization of the weights across particles.
weights = ones(1,number_of_particles)/number_of_particles ;
StateVectors = bsxfun(@plus,StateVectorVarianceSquareRoot*randn(state_variance_rank,number_of_particles),StateVectorMean);
if pruning
    StateVectors_ = StateVectors;
end

% Loop over observations
for t=1:sample_size
    yhat = bsxfun(@minus,StateVectors,state_variables_steady_state);
    epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,number_of_particles);
    if pruning
        yhat_ = bsxfun(@minus,StateVectors_,state_variables_steady_state);
        [tmp, tmp_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,steadystate,DynareOptions.threads.local_state_space_iteration_2);
    else
        tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,DynareOptions.threads.local_state_space_iteration_2);
    end
    PredictedObservedMean = tmp(mf1,:)*transpose(weights);
    PredictionError = bsxfun(@minus,Y(:,t),tmp(mf1,:));
    dPredictedObservedMean = bsxfun(@minus,tmp(mf1,:),PredictedObservedMean);
    PredictedObservedVariance = bsxfun(@times,dPredictedObservedMean,weights)*dPredictedObservedMean' + H;
    lnw = -.5*(const_lik+log(det(PredictedObservedVariance))+sum(PredictionError.*(PredictedObservedVariance\PredictionError),1));
    dfac = max(lnw);
    wtilde = weights.*exp(lnw-dfac);
    lik(t) = log(sum(wtilde))+dfac;
    weights = wtilde/sum(wtilde);
    if (strcmp(DynareOptions.particle.resampling.status,'generic') && neff(weights)<DynareOptions.particle.resampling.neff_threshold*sample_size ) || strcmp(DynareOptions.particle.resampling.status,'systematic')
        idx = resample(weights,DynareOptions.particle.resampling.method1,DynareOptions.particle.resampling.method2);
        StateVectors = tmp(mf0,idx);
        if pruning
            StateVectors_ = tmp_(mf0,idx);
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif strcmp(DynareOptions.particle.resampling.status,'smoothed')
        StateVectors = multivariate_smooth_resampling(weights',tmp(mf0,:)',number_of_particles,DynareOptions.particle.resampling.number_of_partitions)';
        if pruning
            StateVectors_ = multivariate_smooth_resampling(weights',tmp_(mf0,:)',number_of_particles,DynareOptions.particle.resampling.number_of_partitions)';
        end
        weights = ones(1,number_of_particles)/number_of_particles;
    elseif strcmp(DynareOptions.particle.resampling.status,'none')
        StateVectors = tmp(mf0,:);
        if pruning
            StateVectors_ = tmp_(mf0,:);
        end
    end
end

LIK = -sum(lik(start:end));



function n = neff(w)
    n = dot(w,w);