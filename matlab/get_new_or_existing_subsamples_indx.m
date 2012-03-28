function subsamples_indx = get_new_or_existing_subsamples_indx(name1, name2)
% function subsamples_indx = get_new_or_existing_subsamples_indx(name1, name2)
%
% Returns the new estimation_info.subsamples_index index
% for the name1 & name2 pair
%
% INPUTS
%    name1              [string]  the variable for which the subsample was declared
%    name2              [string]  used only in case of corr(name1,name2).subsamples()
%
% OUTPUTS
%    subsamples_indx    [integer] new index in the estimation_info.subsamples
%                                 structure associated with the name pair
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2012 Dynare Team
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

global estimation_info

if isempty(estimation_info.subsamples_index)
    subsamples_indx = 1;
    return;
end

if isempty(name2) % parameter or std() statement
    subsamples_indx = find(strcmp(name1, estimation_info.subsamples_index) == 1);
else % corr statement
    subsamples_indx = find(strcmp([name1 '_' name2], estimation_info.subsamples_index) == 1);
    if isempty(subsamples_indx)
        subsamples_indx = find(strcmp([name2 '_' name1], estimation_info.subsamples_index) == 1);
    end
end

if isempty(subsamples_indx)
    subsamples_indx = size(estimation_info.subsamples_index, 2) + 1;
end

if size(subsamples_indx,2) > 1
    error(['Error: ' name1 ' ' name2 'found more than once in estimation_info.subsamples_index']);
end
