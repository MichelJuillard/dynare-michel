function [LIK,lik] = monte_carlo_gaussian_particle_filter(reduced_form_model,Y,start,number_of_particles)
% hparam,y,nbchocetat,nbchocmesure,smol_prec,nb_part,g,m,choix
% Evaluates the likelihood of a non linear model assuming that the particles are normally distributed. 
%
% INPUTS
%    reduced_form_model     [structure] Matlab's structure describing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    number_of_particles    [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2009-2010 Dynare Team
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

global M_ bayestopt_ oo_
persistent init_flag
persistent restrict_variables_idx observed_variables_idx state_variables_idx mf0 mf1
persistent sample_size number_of_state_variables number_of_observed_variables number_of_structural_innovations

% Set defaults.
if (nargin<4) || (nargin==4 && isempty(number_of_particles))
    number_of_particles = 10 ;
end
if nargin==2 || isempty(start)
   start = 1; 
end

dr = reduced_form_model.state.dr;% Decision rules and transition equations.
Q  = reduced_form_model.state.Q;% Covariance matrix of the structural innovations.
H  = reduced_form_model.measurement.H;% Covariance matrix of the measurement errors.

% Set persistent variables.
if isempty(init_flag)
    mf0 = bayestopt_.mf0;
    mf1 = bayestopt_.mf1;
    restrict_variables_idx  = bayestopt_.restrict_var_list;
    observed_variables_idx  = restrict_variables_idx(mf1);
    state_variables_idx     = restrict_variables_idx(mf0);
    sample_size = size(Y,2);
    number_of_state_variables = length(mf0);
    number_of_observed_variables = length(mf1);
    number_of_structural_innovations = length(Q); 
    init_flag = 1;
end

% Set local state space model (second order approximation).
ghx = dr.ghx(restrict_variables_idx,:);
ghu = dr.ghu(restrict_variables_idx,:);
half_ghxx = .5*dr.ghxx(restrict_variables_idx,:);
half_ghuu = .5*dr.ghuu(restrict_variables_idx,:);
ghxu = dr.ghxu(restrict_variables_idx,:);
steadystate = dr.ys(dr.order_var(restrict_variables_idx));
constant = steadystate + .5*dr.ghs2(restrict_variables_idx);
state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));

StateVectorMean = state_variables_steady_state;
StateVectorVariance = lyapunov_symm(ghx(mf0,:),ghu(mf0,:)*Q*ghu(mf0,:)',1e-12,1e-12);
StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot,2);

Q_lower_triangular_cholesky = chol(Q)';

% Set seed for randn().
seed  = [ 362436069 ; 521288629 ];
randn('state',seed);
 
const_lik = log(2*pi)*number_of_observed_variables; 
lik  = NaN(sample_size,1);

for t=1:sample_size
    PredictedStateMean = zeros(number_of_state_variables,1);
    PredictedObservedMean = zeros(number_of_observed_variables,1);
    PredictedStateVariance = zeros(number_of_state_variables);
    PredictedObservedVariance = zeros(number_of_observed_variables);
    PredictedStateAndObservedCovariance = zeros(number_of_state_variables,number_of_observed_variables);
    for i=1:number_of_particles
        StateVector = StateVectorMean + StateVectorVarianceSquareRoot*randn(state_variance_rank,1);
        yhat = StateVector-state_variables_steady_state;
        epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,1);
        tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,half_ghxx,half_ghuu,ghxu);
        PredictedStateMean = PredictedStateMean + (tmp(mf0))/number_of_particles ;
        PredictedObservedMean = PredictedObservedMean + (tmp(mf1))/number_of_particles;
        PredictedStateVariance = PredictedStateVariance + (tmp(mf0)*tmp(mf0)')/(number_of_particles) ;
        PredictedObservedVariance = PredictedObservedVariance + (tmp(mf1)*tmp(mf1)')/(number_of_particles);
        PredictedStateAndObservedCovariance = PredictedStateAndObservedCovariance + (tmp(mf0)*tmp(mf1)')/(number_of_particles) ;
    end
    PredictedObservedVariance = PredictedObservedVariance + H - PredictedObservedMean*(PredictedObservedMean');
    PredictedStateVariance = PredictedStateVariance - PredictedStateMean*(PredictedStateMean') ;
    PredictedStateAndObservedCovariance = PredictedStateAndObservedCovariance - PredictedStateMean*(PredictedObservedMean');
    iPredictedObservedVariance = inv(PredictedObservedVariance);
    prediction_error = Y(:,t) - PredictedObservedMean;
    filter_gain = PredictedStateAndObservedCovariance*iPredictedObservedVariance;
    StateVectorMean = PredictedStateMean + filter_gain*prediction_error;
    StateVectorVariance = PredictedStateVariance - filter_gain*PredictedObservedVariance*filter_gain';
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)';
    state_variance_rank = size(StateVectorVarianceSquareRoot,2);
    lik(t) = -.5*(const_lik + log(det(PredictedObservedVariance)) + prediction_error'*iPredictedObservedVariance*prediction_error);
end

LIK = -sum(lik(start:end));