function mhname = get_name_of_the_last_mh_file(M_)
%function mhname = get_name_of_the_last_mh_file(M_)

% Copyright (C) 2008 Dynare Team
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

    model_name = M_.fname ;
    mcmc_directory = M_.dname ;
    load([ mcmc_directory '/metropolis/' model_name '_mh_history.mat']) ;
    mh_number = record.LastFileNumber ;
    bk_number = record.Nblck ;
    clear('record') ;
    mhname = [ mcmc_directory ...
               '/metropolis/' ...
               model_name ...
               '_mh' ...
               int2str(mh_number) ...
               '_blck' ...
               int2str(bk_number) ...
               '.mat'] ;