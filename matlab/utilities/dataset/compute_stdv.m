function dataset_ = compute_stdv(dataset_) 
% Compute the standard deviation for each observed variable (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{dataset_} =} compute_stdv(@var{dataset_})
%! @anchor{compute_stdv}
%! This function computes the standard deviation of the observed variables (possibly with missing observations).
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
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.stdv} is a 
%! @tex{n\times 1} vector (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}).
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

if dataset_.missing.state
    dataset_.descriptive.stdv = sqrt(nanmean(bsxfun(@power,nandemean(transpose(dataset_.data)),2)));
else
    dataset_.descriptive.stdv = sqrt(mean(bsxfun(@power,demean(transpose(dataset_.data)),2)));
end