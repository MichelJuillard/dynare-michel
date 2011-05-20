function [options_, oo_]=ms_forecast(M_, options_, oo_)
%function ms_forecast()
% calls ms irf mex function
%
% INPUTS
%    M_
%    options_
%    oo_
%
% OUTPUTS
%    options_
%    oo_
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011 Dynare Team
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

disp('Forecasts');
options_ = set_ms_estimation_flags_for_other_mex(options_);
options_ = set_ms_simulation_flags_for_other_mex(options_);
oo_ = set_oo_w_estimation_output(options_, oo_);

opt = {options_.ms.output_file_tag, ...
    'seed', options_.DynareRandomStreams.seed, ...
    'horizon', options_.ms.horizon, ...
    'data',options_.ms.forecast_data_obs ...
    'error_bands', options_.ms.error_bands, ...
    'percentiles', options_.ms.error_band_percentiles, ...
    'thin', options_.ms.thinning_factor };

[err, forecast] = mex_ms_forecast(opt{:},'free_parameters',oo_.ms.maxparams,'shocks_per_parameter', options_.ms.shock_draws);
mexErrCheck('mex_ms_forecast ergodic ', err);
plot_ms_forecast(M_,forecast,'Forecast');

[err, regime_forecast] = mex_ms_forecast(opt{:},'free_parameters',oo_.ms.maxparams,'shocks_per_parameter', options_.ms.shock_draws,'regimes');
mexErrCheck('mex_ms_forecast ergodic regimes', err);
save([M_.fname '/' options_.ms.output_file_tag '_ergodic_forecast.mat'], 'forecast', 'regime_forecast');

if exist(options_.ms.load_mh_file,'file') > 0
    [err, forecast] = mex_ms_forecast(opt{:},'free_parameters',oo_.ms.maxparams,'shocks_per_parameter', options_.ms.shocks_per_parameter, ...
        'simulation_file',options_.ms.load_mh_file,'parameter_uncertainty');
    mexErrCheck('mex_ms_forecast bayesian ', err);
    plot_ms_forecast(M_,forecast,'Forecast w/ Parameter Uncertainty');

    [err, regime_forecast] = mex_ms_forecast(opt{:},'free_parameters',oo_.ms.maxparams,'shocks_per_parameter', options_.ms.shocks_per_parameter, ...
        'simulation_file',options_.ms.load_mh_file,'parameter_uncertainty','regimes');
    mexErrCheck('mex_ms_forecast bayesian regimes ', err);
    save([M_.fname '/' options_.ms.output_file_tag '_bayesian_forecast.mat'], 'forecast', 'regime_forecast');
end
options_ = initialize_ms_sbvar_options(M_, options_);
end
