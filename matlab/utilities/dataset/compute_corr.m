function dataset_ = compute_corr(dataset_) 
% Computes the correlation matrix of the sample (possibly with missing observations).

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
%! @ref{ndim}, @ref{compute_cova}. 
%!    
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.corr} is a 
%! @tex{n\times n} vector (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}).
%!
%! @strong{Remark 2.} If @code{dataset_.descriptive.cova} does not exist, the covariance matrix is computed prior to the 
%! computation of the correlation matrix.
%!    
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2012 Dynare Team
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

if ~isfield(dataset_.descriptive,'cova')
    dataset_ = compute_cova(dataset_);
end
normalization_matrix = diag(1./sqrt(diag(dataset_.descriptive.cova)));
dataset_.descriptive.corr = normalization_matrix*dataset_.descriptive.cova*normalization_matrix;