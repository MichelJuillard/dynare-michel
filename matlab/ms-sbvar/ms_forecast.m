function [options_, oo_]=ms_forecast(M_, options_, oo_)
%function ms_forecast()
% MS-SBVAR Forecast
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
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

disp('MS-SBVAR Forecasts');
options_ = set_file_tags(options_);
[options_, oo_] = set_ms_estimation_file(options_, oo_);
options_ = set_ms_simulation_file(options_);
clean_files_for_second_type_of_mex(M_, options_, 'forecast')
forecastdir = [options_.ms.output_file_tag filesep 'Forecast'];
create_dir(forecastdir);

opt = { ...
    {'file_tag', options_.ms.estimation_file_tag}, ...
    {'seed', options_.DynareRandomStreams.seed}, ...
    {'horizon', options_.ms.horizon}, ...
    {'number_observations', options_.ms.forecast_data_obs}, ...
    {'error_bands', options_.ms.error_bands}, ...
    {'percentiles', options_.ms.error_band_percentiles}, ...
    {'thin', options_.ms.thinning_factor}
    };

if options_.ms.median
    opt = [opt(:)' {{'median'}}];
end

[err, forecast] = mex_ms_forecast([opt(:)', {{'free_parameters',oo_.ms.maxparams}, ...
    {'shocks_per_parameter', options_.ms.shock_draws}}]);
mexErrCheck('mex_ms_forecast ergodic ', err);
plot_ms_forecast(M_,options_,forecast,'Forecast',options_.graph_save_formats,options_.TeX);

[err, regime_forecast] = mex_ms_forecast([opt(:)', {{'free_parameters',oo_.ms.maxparams}, ...
    {'shocks_per_parameter', options_.ms.shock_draws}, {'regimes'}}]);
mexErrCheck('mex_ms_forecast ergodic regimes', err);
save([forecastdir filesep 'ergodic_forecast.mat'], 'forecast', 'regime_forecast');

if exist(options_.ms.mh_file,'file') > 0
    [err, forecast] = mex_ms_forecast([opt(:)', {{'free_parameters',oo_.ms.maxparams}, ...
        {'shocks_per_parameter', options_.ms.shocks_per_parameter}, ...
        {'simulation_file', options_.ms.mh_file}, {'parameter_uncertainty'}}]);
    mexErrCheck('mex_ms_forecast bayesian ', err);
    plot_ms_forecast(M_,options_,forecast,'Forecast w/ Parameter Uncertainty',options_.graph_save_formats,options_.TeX);

    [err, regime_forecast] = mex_ms_forecast([opt(:)', {{'free_parameters',oo_.ms.maxparams}, ...
        {'shocks_per_parameter', options_.ms.shocks_per_parameter}, ...
        {'simulation_file', options_.ms.mh_file}, {'parameter_uncertainty','regimes'}}]);
    mexErrCheck('mex_ms_forecast bayesian regimes ', err);
    save([forecastdir filesep 'bayesian_forecast.mat'], 'forecast', 'regime_forecast');
end
end
