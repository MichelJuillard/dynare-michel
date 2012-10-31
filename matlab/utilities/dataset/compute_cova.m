function dataset_ = compute_cova(dataset_) 
% Computes the covariance matrix of the sample (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{dataset_} =} compute_corr(@var{dataset_})
%! @anchor{compute_corr}
%! This function computes covariance matrix of the sample (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item dataset_
%! Dynare structure describing the dataset, built by @ref{initialize_dataset}
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item dataset_
%! Dynare structure describing the dataset, built by @ref{initialize_dataset}
%! @end table
%! 
%! @strong{This function is called by:} 
%! @ref{descriptive_statistics}.
%! 
%! @strong{This function calls:}
%! @ref{ndim}, @ref{demean}, @ref{nandemean}.
%!    
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.cova} is a 
%! @tex{n\times n} vector (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}).
%!    
%! @end deftypefn
%@eod:
    
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

dataset_.descriptive.cova = zeros(dataset_.nvobs);

data = transpose(dataset_.data);

for i=1:dataset_.info.nvobs
    for j=i:dataset_.info.nvobs
        if dataset_.missing.state
            dataset_.descriptive.cova(i,j) = nanmean(nandemean(data(:,i)).*nandemean(data(:,j)));
        else
            dataset_.descriptive.cova(i,j) = mean(demean(data(:,i)).*demean(data(:,j)));
        end
        if j>i
            dataset_.descriptive.cova(j,i) = dataset_.descriptive.cova(i,j);
        end
    end
end