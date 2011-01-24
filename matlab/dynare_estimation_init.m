function [data,rawdata]=dynare_estimation_init(var_list_, igsa)

% function dynare_estimation_init(var_list_)
% preforms initialization tasks before estimation or
% global sensitivity analysis
%  
% INPUTS
%   var_list_:  selected endogenous variables vector
%  
% OUTPUTS
%   data:    data after required transformation
%   rawdat:  data as in the data file
%
% SPECIAL REQUIREMENTS
%   none

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

global M_ options_ oo_ estim_params_ 
global bayestopt_

if nargin<2 | isempty(igsa),
    igsa=0;
end

options_.varlist = var_list_;
options_.lgyidx2varobs = zeros(size(M_.endo_names,1),1);
for i = 1:size(M_.endo_names,1)
    tmp = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
    if ~isempty(tmp)
        options_.lgyidx2varobs(i,1) = tmp;
    end
end

if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    options_.dsge_var = 1;
end

if options_.order > 1
    options_.order = 1;
end

if options_.prefilter == 1
    options_.noconstant = 1;
end

if options_.filtered_vars ~= 0 & options_.filter_step_ahead == 0
    options_.filter_step_ahead = 1;
end
if options_.filter_step_ahead ~= 0
    options_.nk = max(options_.filter_step_ahead);
else
    options_.nk = 0;
end

%% Add something to the parser ++>
% The user should be able to choose another name
% for the directory...
M_.dname = M_.fname; 

pnames          = ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
n_varobs        = size(options_.varobs,1);

if ~isempty(estim_params_)
    [xparam1,estim_params_,bayestopt_,lb,ub, M_] = set_prior(estim_params_, M_, options_);

    if any(bayestopt_.pshape > 0)
        if options_.mode_compute
            plot_priors
        end
    else
        options_.mh_replic = 0;
    end

    % set prior bounds and check initial value of the parameters
    bounds = prior_bounds(bayestopt_);
    bounds(:,1)=max(bounds(:,1),lb);
    bounds(:,2)=min(bounds(:,2),ub);

    if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
        find(xparam1 < bounds(:,1))
        find(xparam1 > bounds(:,2))
        error('Initial parameter values are outside parameter bounds')
    end
    lb = bounds(:,1);
    ub = bounds(:,2);
    bayestopt_.lb = lb;
    bayestopt_.ub = ub;
else
    xparam1 = [];
    bayestopt_.lb = [];
    bayestopt_.ub = [];
    bayestopt_.jscale = [];
    bayestopt_.pshape = [];
    bayestopt_.p1 = [];
    bayestopt_.p2 = [];
    bayestopt_.p3 = [];
    bayestopt_.p4 = [];
    estim_params_.nvx = 0;
    estim_params_.nvn = 0;
    estim_params_.ncx = 0;
    estim_params_.ncn = 0;
    estim_params_.np = 0;
end
nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np ;
nx  = nvx+nvn+ncx+ncn+np;

if ~isfield(options_,'trend_coeffs')
    bayestopt_.with_trend = 0;
else
    bayestopt_.with_trend = 1;
    bayestopt_.trend_coeff = {};
    trend_coeffs = options_.trend_coeffs;
    nt = length(trend_coeffs);
    for i=1:n_varobs
        if i > length(trend_coeffs)
            bayestopt_.trend_coeff{i} = '0';
        else
            bayestopt_.trend_coeff{i} = trend_coeffs{i};
        end
    end
end

bayestopt_.penalty = 1e8;       % penalty 

dr = set_state_space([],M_);
nstatic = dr.nstatic;
npred = dr.npred;
nspred = dr.nspred;

if isempty(options_.varobs)
    error('ESTIMATION: VAROBS is missing')
end

%% Setting resticted state space (observed + predetermined variables)

k = [];
k1 = [];
for i=1:n_varobs
    k = [k strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:),'exact')];
    k1 = [k1 strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
end
% union of observed and state variables
k2 = union(k',[dr.nstatic+1:dr.nstatic+dr.npred]', 'rows');

% set restrict_state to postion of observed + state variables
% in expanded state vector
bayestopt_.restrict_var_list = k2;
% set mf0 to positions of state variables in restricted state vector
% for likelihood computation.
[junk,bayestopt_.mf0] = ismember([dr.nstatic+1:dr.nstatic+dr.npred]',k2);
% set mf1 to positions of observed variables in restricted state vector
% for likelihood computation.
[junk,bayestopt_.mf1] = ismember(k,k2);
% set mf2 to positions of observed variables in expanded state vector
% for filtering and smoothing
bayestopt_.mf2  = k;
bayestopt_.mfys = k1;

[junk,ic] = intersect(k2,nstatic+(1:npred)');
bayestopt_.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];
aux = dr.transition_auxiliary_variables;
aux(:,2) = aux(:,2) + sum(k2 <= nstatic);
k = find(aux(:,2) > npred);
aux(k,2) = aux(k,2) + sum(k2 > nstatic+npred);
bayestopt_.restrict_aux = aux;


%% Initialization with unit-root variables
if ~isempty(options_.unit_root_vars)
    n_ur = size(options_.unit_root_vars,1);
    i_ur = zeros(n_ur,1);
    for i=1:n_ur
        i1 = strmatch(deblank(options_.unit_root_vars(i,:)),M_.endo_names(dr.order_var,:),'exact');
        if isempty(i1)
            error('Undeclared variable in unit_root_vars statement')
        end
        i_ur(i) = i1;
    end
    bayestopt_.var_list_stationary = setdiff((1:M_.endo_nbr)',i_ur);
    [junk,bayestopt_.restrict_var_list_nonstationary] = ...
        intersect(bayestopt_.restrict_var_list,i_ur);
    bayestopt_.restrict_var_list_stationary = ...
        setdiff((1:length(bayestopt_.restrict_var_list))', ...
                bayestopt_.restrict_var_list_nonstationary);
    if M_.maximum_lag > 1
        l1 = flipud([cumsum(M_.lead_lag_incidence(1:M_.maximum_lag-1,dr.order_var),1);ones(1,M_.endo_nbr)]);
        l2 = l1(:,bayestopt_.restrict_var_list);
        il2 = find(l2' > 0);
        l2(il2) = (1:length(il2))';
        bayestopt_.restrict_var_list_stationary = ...
            nonzeros(l2(:,bayestopt_.restrict_var_list_stationary)); 
        bayestopt_.restrict_var_list_nonstationary = ...
            nonzeros(l2(:,bayestopt_.restrict_var_list_nonstationary)); 
    end
    options_.lik_init = 3;
end % if ~isempty(options_.unit_root_vars)

if isempty(options_.datafile),
    if igsa,
        data=[];
        rawdata=[];
        return,
    else
        error('ESTIMATION: datafile option is missing'),
    end
end

%% If jscale isn't specified for an estimated parameter, use
%% global option options_.jscale, set to 0.2, by default
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

%% Read and demean data 
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);

options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1 & ~options_.logdata
    rawdata = log(rawdata);
end
if options_.prefilter == 1
    bayestopt_.mean_varobs = mean(rawdata,1)';
    data = transpose(rawdata-repmat(bayestopt_.mean_varobs',gend,1));
else
    data = transpose(rawdata);
end

if ~isreal(rawdata)
    error(['There are complex values in the data. Probably  a wrong' ...
           ' transformation'])
end

