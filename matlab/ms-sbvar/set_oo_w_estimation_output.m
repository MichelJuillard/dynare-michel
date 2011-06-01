function oo_=set_oo_w_estimation_output(options_, oo_)
%function set_oo_w_estimation_output()
% places estimation output in oo_ structure
%
% INPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
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


oo_.ms.maxparams = load(options_.ms.free_param_file);
oo_.ms.maxparams = oo_.ms.maxparams(3:end)';

if ~isfield(oo_.ms, 'A0') || ~isfield(oo_.ms, 'Aplus') ...
        || ~isfield(oo_.ms, 'Zeta') || ~isfield(oo_.ms, 'Q')
    [err, oo_.ms.A0, oo_.ms.Aplus, oo_.ms.Zeta, oo_.ms.Q] = ...
        mex_ms_convert_free_parameters(options_.ms.output_file_tag, oo_.ms.maxparams);
    mexErrCheck('mex_ms_convert_free_parameters', err);
end
end
