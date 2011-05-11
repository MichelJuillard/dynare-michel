function clean_ms_files(output_file_tag)
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

delete_if_exist(['est_aux_' output_file_tag '.out']);
delete_if_exist(['est_csminwel_' output_file_tag '.out']);
delete_if_exist(['est_final_' output_file_tag '.out']);
delete_if_exist(['est_flat_header_' output_file_tag '.out']);
delete_if_exist(['est_flat_' output_file_tag '.out']);
delete_if_exist(['est_free_' output_file_tag '.out']);
delete_if_exist(['est_intermediate_' output_file_tag '.out']);
delete_if_exist('g1.mat');
delete_if_exist('H.dat');
delete_if_exist(['init_' output_file_tag '.dat']);
delete_if_exist(['matlab_' output_file_tag '.prn']);
delete_if_exist(['mdd_t3_' output_file_tag '.out']);
delete_if_exist(['proposal_t3_' output_file_tag '.out']);
delete_if_exist(['simulation_info_' output_file_tag '.out']);
delete_if_exist(['simulation_' output_file_tag '.out']);
delete_if_exist([output_file_tag '_markov_file.dat']);
end

function delete_if_exist(fname)
if exist(fname,'file') == 2
    delete(fname);
end
end
