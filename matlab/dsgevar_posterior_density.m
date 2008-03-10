function bvar = dsgevar_posterior_density(deep)
% This function characterizes the posterior distribution of a bvar with 
% a dsge prior (as in Del Negro and Schorfheide 2003) for a given value 
% of the deep parameters (structural parameters + the size of the 
% shocks + dsge_prior_weight).
%    
% INPUTS
%   deep:      [double] a vector with the deep parameters.
%  
% OUTPUTS
%   bvar:      a matlab structure with prior and posterior densities. 
%  
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none
%  
%  
% part of DYNARE, copyright Dynare Team (1996-2007)
% Gnu Public License.
global options_ M_

gend = options_.nobs;
dsge_prior_weight = M_.params(strmatch('dsge_prior_weight',M_.param_names));
DSGE_PRIOR_WEIGHT = floor(gend*(1+dsge_prior_weight));

bvar.NumberOfLags = options_.varlag;
bvar.NumberOfVariables = size(options_.varobs,1);
bvar.Constant = 'no';
bvar.NumberOfEstimatedParameters = bvar.NumberOfLags*bvar.NumberOfVariables;   
if 
    bvar.Constant = 'yes';
    bvar.NumberOfEstimatedParameters = bvar.NumberOfEstimatedParameters + ...
        bvar.NumberOfVariables;        
end

[fval,cost_flag,info,PHI,SIGMAu,iXX,prior] =  DsgeVarLikelihood(deep',gend);

% Conditionnal posterior density of the lagged matrices (given Sigma) ->
% Matric-variate normal distribution.
bvar.LaggedMatricesConditionalOnSigma.posterior.density = 'matric-variate normal';
bvar.LaggedMatricesConditionalOnSigma.posterior.arg1 = PHI;
bvar.LaggedMatricesConditionalOnSigma.posterior.arg2 = 'Sigma';
bvar.LaggedMatricesConditionalOnSigma.posterior.arg3 = iXX;

% Marginal posterior density of the covariance matrix -> Inverted Wishart. 
bvar.Sigma.posterior.density = 'inverse wishart';
bvar.Sigma.posterior.arg1 = SIGMAu*DSGE_PRIOR_WEIGHT;
bvar.Sigma.posterior.arg2 = DSGE_PRIOR_WEIGHT-bvar.NumberOfEstimatedParameters;

% Marginal posterior density of the lagged matrices -> Generalized 
% Student distribution (See appendix B.5 in Zellner (1971)).
bvar.LaggedMatrices.posterior.density = 'matric-variate student';
bvar.LaggedMatrices.posterior.arg1 = inv(iXX);%P
bvar.LaggedMatrices.posterior.arg2 = SIGMAu*DSGE_PRIOR_WEIGHT;%Q
bvar.LaggedMatrices.posterior.arg3 = PHI;%M (posterior mean)
bvar.LaggedMatrices.posterior.arg4 = DSGE_PRIOR_WEIGHT;%(sample size)



% Conditionnal posterior density of the lagged matrices (given Sigma) ->
% Matric-variate normal distribution.
bvar.LaggedMatricesConditionalOnSigma.prior.density = 'matric-variate normal';
bvar.LaggedMatricesConditionalOnSigma.prior.arg1 = prior.PHIstar;
bvar.LaggedMatricesConditionalOnSigma.prior.arg2 = 'Sigma';
bvar.LaggedMatricesConditionalOnSigma.prior.arg3 = prior.iGXX;

% Marginal posterior density of the covariance matrix -> Inverted Wishart. 
bvar.Sigma.prior.density = 'inverse wishart';
bvar.Sigma.prior.arg1 = prior.SIGMAstar*prior.ArtificialSampleSize;
bvar.Sigma.prior.arg2 = prior.DF;

% Marginal posterior density of the lagged matrices -> Generalized 
% Student distribution (See appendix B.5 in Zellner (1971)).
bvar.LaggedMatrices.prior.density = 'matric-variate student';
bvar.LaggedMatrices.prior.arg1 = inv(iGXX);%P
bvar.LaggedMatrices.prior.arg2 = prior.SIGMAstar*prior.ArtificialSampleSize;%Q
bvar.LaggedMatrices.prior.arg3 = prior.PHIstar;%M (posterior mean)
bvar.LaggedMatrices.prior.arg4 = prior.ArtificialSampleSize;%(sample size)