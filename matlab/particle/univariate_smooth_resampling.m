function new_particles = univariate_smooth_resmp(weights,particles,number)
% Smooth Resamples particles.

%@info:
%! @deftypefn {Function File} {@var{indx} =} resample (@var{weights}, @var{method})
%! @anchor{particle/resample}
%! @sp 1
%! Resamples particles.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item weights
%! n*1 vector of doubles, particles' weights.
%! @item method
%! string equal to 'residual' or 'traditional'.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item indx
%! n*1 vector of intergers, indices.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{particle/sequantial_importance_particle_filter}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{residual_resampling}, @ref{traditional_resampling}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2011, 2012 Dynare Team
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

% AUTHOR(S) frederic DOT karame AT univ DASH evry DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

M = length(particles) ;
lambda_tilde = [ (.5*(2*weights(1)+weights(2))) ;       
                 (.5*(weights(2:M-1)+weights(3:M))) ;
                 (.5*(weights(M-1)+2*weights(M))) ] ;
lambda_bar = cumsum(lambda_tilde) ;
lambda_bar = lambda_bar(1:M-1) ;
u = rand(1,1) ;
new_particles = zeros(number,1) ;
i = 1 ;
j = 1 ;
while i<=number 
    u_j = ( i-1 + u)/number ;
    while u_j>lambda_bar(j) 
        j = j+1 ;
        if j==M 
            j = M-1 ;
            break ;
        end 
    end  
    u_star = (u_j - (lambda_bar(j)-lambda_tilde(j)))./lambda_tilde(j) ;
    new_particles(i) = (particles(j+1) - particles(j))*u_star + particles(j) ;
    i = i+1 ;
end 

