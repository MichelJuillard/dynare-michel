function new_particles = multivariate_smooth_resampling(weights,particles,number_of_new_particles,number_of_partitions)
% Smooth Resampling of the  particles.

%@info:
%! @deftypefn {Function File} {@var{new_particles} =} multivariate_smooth_resampling (@var{weights}, @var{particles}, @var{number_of_new_particles}, @var{number_of_partitions})
%! @anchor{particle/multivariate_smooth_resampling}
%! @sp 1
%! Smooth Resampling of the  particles (multivariate version).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item weights
%! n*1 vector of doubles, particles' weights.
%! @item particles
%! n*1 vector of doubles, particles.
%! @item number_of_new_particles
%! Integer scalar.
%! @item number_of_partitions
%! Integer scalar.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item indx
%! number_of_new_particles*1 vector of doubles, new particles.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{particle/sequantial_importance_particle_filter}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{particle/univariate_smooth_resampling}
%! @sp 2
%! @end deftypefn
%@eod:

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

% AUTHOR(S) frederic DOT karame AT univ DASH lemans DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

number_of_particles = length(weights);
number_of_states = size(particles,2);
number = number_of_particles/number_of_partitions ;
tout = sortrows([particles weights],1) ;
particles = tout(:,1:number_of_states) ;
weights = tout(:,1+number_of_states) ;
if number_of_partitions>1 
    cum_weights = cumsum(weights) ;
    indx = 1:number_of_particles ;
    for i=1:number_of_partitions
        if i==number_of_partitions
          tmp = bsxfun(@gt,cum_weights,(i-1)/number_of_partitions) ;
          kp = indx( tmp ) ;
        else
          tmp = bsxfun(@and,bsxfun(@gt,cum_weights,(i-1)/number_of_partitions),bsxfun(@lt,cum_weights,i/number_of_partitions)) ;
          kp = indx( tmp ) ;
        end
        if numel(kp)>2
            Np = length(kp) ;
            wtilde = [ ( number_of_partitions*( cum_weights(kp(1)) - (i-1)/number_of_partitions) ) ;
                       ( number_of_partitions*weights(kp(2:Np-1)) ) ;
                       ( number_of_partitions*(i/number_of_partitions - cum_weights(kp(Np)-1) ) ) ] ;
        elseif numel(kp)==2
            Np = length(kp) ;
            wtilde = [ ( number_of_partitions*( cum_weights(kp(1)) - (i-1)/number_of_partitions) ) ;
                       ( number_of_partitions*(i/number_of_partitions - cum_weights(kp(Np)-1) ) ) ] ;        
        elseif numel(kp)==1 
            new_particles = ones(number,1).*particles(kp,:) ;    
            return ;
        else
            % probleme
        end
        test = sum(wtilde) ;
        disp(test)
        new_particles = zeros(number_of_particles,number_of_states) ;
        new_particles_j = zeros(number,number_of_states) ;
        for j=1:number_of_states
          particles_j = particles(kp,j) ;
          if j>1
            tout = sortrows( [ particles_j wtilde],1) ;
            particles_j = tout(:,1) ;
            wtilde = tout(:,2) ;
          end
          new_particles_j(:,j) = univariate_smooth_resampling(wtilde,particles_j,number) ;
        end
        new_particles((i-1)*number+1:i*number,:) = new_particles_j;
    end
else
    new_particles = zeros(number,number_of_states) ;
    for j=1:number_of_states
      if j>1
        tout = sortrows( [ particles(:,j) weights],1) ;
        particles_j = tout(:,1) ;
        weights_j = tout(:,2) ;
      else
        particles_j = particles(:,j) ;
        weights_j = weights ;
      end
      new_particles(:,j) = univariate_smooth_resampling(weights_j,particles_j,number) ;
    end
end