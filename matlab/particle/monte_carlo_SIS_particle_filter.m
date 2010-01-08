function [LIK,lik] = monte_carlo_SIS_particle_filter(reduced_form_model,Y,start,number_of_particles)
% hparam,y,nbchocetat,nbchocmesure,smol_prec,nb_part,g,m,choix
% Evaluates the likelihood of a nonlinear model with a particle filter without systematic resampling. 
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
lik = NaN(sample_size,1);
nb_obs_resamp = 0 ;
w = ones(number_of_particles,1) ;
for t=1:sample_size
    PredictedState = zeros(number_of_particles,number_of_state_variables);
    PredictionError = zeros(number_of_particles,number_of_observed_variables);
    %PredictedStateMean = zeros(number_of_state_variables,1);
    PredictedObservedMean = zeros(number_of_observed_variables,1);
    %PredictedStateVariance = zeros(number_of_state_variables,number_of_state_variables);
    PredictedObservedVariance = zeros(number_of_observed_variables,number_of_observed_variables);
    %PredictedStateAndObservedCovariance = zeros(number_of_state_variables,number_of_observed_variables);
    for i=1:number_of_particles
        if t==1
          StateVector = StateVectorMean + StateVectorVarianceSquareRoot*randn(state_variance_rank,1);
        else 
          StateVector = StateUpdated(i,:)' ;  
        end 
        yhat = StateVector-state_variables_steady_state;
        epsilon = Q_lower_triangular_cholesky*randn(number_of_structural_innovations,1);
        tmp = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,half_ghxx,half_ghuu,ghxu);
        % stockage des particules et des erreurs de prévisions
        PredictedState(i,:) = tmp(mf0)' ;
        PredictionError(i,:) = (Y(:,t) - tmp(mf1))' ; 
        %PredictedStateMean_old = PredictedStateMean;
        PredictedObservedMean_old = PredictedObservedMean;
        %PredictedStateMean = PredictedStateMean + (tmp(mf0)-PredictedStateMean)/i;
        PredictedObservedMean = PredictedObservedMean + (tmp(mf1)-PredictedObservedMean)/i;
        %psm = PredictedStateMean*PredictedStateMean';
        pom = PredictedObservedMean*PredictedObservedMean';
        %pcm = PredictedStateMean*PredictedObservedMean';
        %PredictedStateVariance = PredictedStateVariance ...
        %    + ( (tmp(mf0)*tmp(mf0)'-psm-PredictedStateVariance)+(i-1)*(PredictedStateMean_old*PredictedStateMean_old'-psm) )/i;
        PredictedObservedVariance = PredictedObservedVariance ...
            + ( (tmp(mf1)*tmp(mf1)'-pom-PredictedObservedVariance)+(i-1)*(PredictedObservedMean_old*PredictedObservedMean_old'-pom) )/i;
        %PredictedStateAndObservedCovariance = PredictedStateAndObservedCovariance ...
        %    + ( (tmp(mf0)*tmp(mf1)'-pcm-PredictedStateAndObservedCovariance)+(i-1)*(PredictedStateMean_old*PredictedObservedMean_old'-pcm) )/i;
    end
    PredictedObservedVariance = PredictedObservedVariance + H;
    iPredictedObservedVariance = inv(PredictedObservedVariance);
    lnw = -0.5*(const_lik + log(det(PredictedObservedVariance)) + sum((PredictionError*iPredictedObservedVariance).*PredictionError,2)) ;
    %bidouille numérique Schorfheide
    dfac = max(lnw);
    wtilde = w.*exp(lnw - dfac) ;
    % vraisemblance de l'observation
    lik(t) = log(mean(wtilde)) + dfac ;
    %clear (PredictionError) ;  
    %clear (lnw) ;
    % calcul des poids 
    w = wtilde/sum(wtilde) ;
    %clear (wtilde) ;
    %update 
    Neff = 1/sum(w.*w) ; 
    if Neff>number_of_particles                        %no resampling
        StateUpdated = PredictedState ; 
        %clear (PredictedState) ;
        w = number_of_particles*w ;
    else                                                %resampling
        nb_obs_resamp = nb_obs_resamp+1 ;

        %kill the smallest particles before resampling :! facultatif ? 
        %to_kill = [w PredictedState] ; 
        %to_kill = delif(to_kill,w<(1/number_of_particules)*1E-12);%%
        %[n,m] = size(to_kill) ;
        %w = to_kill(:,1) ;
        %PredictedState = to_kill(:,2:m) ;
        %clear (to_kill) ;
        %if number_of_particles neq n 
        %  'Elimination de '; number_of_particles - n ; ' particules à l''observation ';t ;
        %end 
        %fin de kill
        %remise à l'échelle des poids sur les particules restantes 
        %w = cumsum( w/sum(w) );
        %Rééchantillonage systématique 
        %rnduvec = ( (1:number_of_particles)-1+rand )/number_of_particles ;
        %selind = (number_of_particles - sum( w > rnduvec ) + 1)'; % problème de mémoire car w .> rnduvec' très grande !
        %clear (rnduvec) ;
        %StateUpdated = PredictedState(selind,:) ;
        %clear (selind) ;
        % initialize
        selind = zeros(number_of_particles,1);
        % construct CDF
        c = cumsum(w);
        % draw a starting point
        rnduvec = ( (1:number_of_particles)-1+rand)/number_of_particles ;
        % start at the bottom of the CDF
        j=1;
        for i=1:number_of_particles
            % move along the CDF
            while (rnduvec(i)>c(j))
                j=j+1;
            end
            % assign index
            selind(i) = j;
        end
        StateUpdated = PredictedState(selind,:);
        w = ones(number_of_particles,1) ;
    end     
end
LIK = -sum(lik(start:end));