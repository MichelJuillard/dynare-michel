function dataset_ = compute_acov(dataset_) 
% Computes the (multivariate) auto-covariance function of the sample (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{dataset_} =} compute_corr(@var{dataset_},@var{nlag})
%! @anchor{compute_acov}
%! This function computes the (multivariate) auto-covariance function of the sample (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item dataset_
%! Dynare structure describing the dataset, built by @ref{initialize_dataset}
%! @item nlag
%! Integer scalar. The maximum number of lags of the autocovariance function.    
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
%! @ref{ndim}, @ref{compute_cova}, @ref{demean}, @ref{nandemean}.
%!    
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.acov} is a 
%! @tex{n\times n\times p} array (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}, 
%! and @tex{n} is the maximum number of lags given by the second input @code{nlag}).
%!    
%! @strong{Remark 2.} If @code{dataset_.descriptive.cova} does not exist, the covariance matrix is computed prior to the 
%! computation of the auto-covariance function.
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

if ~isfield(dataset_.descriptive,'cova')
    dataset_ = compute_cova(dataset_);
end
dataset_.descriptive.acov = zeros(dataset_.nvobs,dataset_.nvobs,nlag);

data = transpose(dataset_.data);

for lag=1:nlag
    for i=1:dataset_.info.nvobs
        for j=1:dataset_.info.nvobs
            if dataset_.missing.state
                dataset_.descriptive.acov(i,j,lag) = nanmean(nandemean(data(lag+1:end,i)).*nandemean(data(1:end-lag,j)));
            else
                dataset_.descriptive.acov(i,j,lag) = mean(demean(data(lag+1:end,i)).*demean(data(1:end-lag,j)));
            end
        end
    end
end