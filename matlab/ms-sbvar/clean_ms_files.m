function clean_ms_files(mod_name)
% function clean_ms_files()
% removes MS intermediary files
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2007-2011 Dynare Team
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

    delete_if_exist(['./draws_' mod_name '.dat'])
    delete_if_exist(['est_aux_' mod_name '.dat'])
    delete_if_exist(['est_csminwel_' mod_name '.dat'])
    delete_if_exist(['est_final_' mod_name '.dat'])
    delete_if_exist(['est_flat_header_' mod_name '.dat'])
    delete_if_exist(['est_flat_' mod_name '.dat'])
    delete_if_exist(['est_intermediate_' mod_name '.dat'])
    delete_if_exist(['header_' mod_name '.dat'])
    delete_if_exist(['init_' mod_name '.dat'])
    delete_if_exist('markov_file.dat')
    delete_if_exist(['matlab_' mod_name '.prn'])
    delete_if_exist(['mhm_draws_' mod_name '.dat'])
    delete_if_exist(['mhm_input_' mod_name '.dat'])
    delete_if_exist(['mhm_intermediate_draws_' mod_name '.dat'])
    delete_if_exist(['mhm_intermediate_' mod_name '.dat'])
    delete_if_exist(['mhm_regime_counts_' mod_name '.dat'])
    delete_if_exist(['probabilities_' mod_name '.dat'])
    delete_if_exist(['truncatedpower_md_posterior_' mod_name '.dat'])
    delete_if_exist(['truncatedpower_md_proposal_' mod_name '.dat'])
    delete_if_exist(['truncatedpower_md_' mod_name '.dat'])

function delete_if_exist(fname)
    if exist(fname) == 2
        delete(fname)
    end
    