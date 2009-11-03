function [LIK,lik] = gaussian_particle_filter(reduced_form_model,Y,start,mf,number_of_particles,grid_size)
% hparam,y,nbchocetat,nbchocmesure,smol_prec,nb_part,g,m,choix
% Evaluates the likelihood of a non linear model assuming that the particles are normally distributed. 
%
% INPUTS
%    reduced_form_model     [structure] Matlab's structure desvcribing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    number_of_particles    [integer]   scalar.
%    grid_size              [integer]   scalar, size of the smoliak grid.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2009 Dynare Team
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

global M_

dr = reduced_form_model.state.dr;
Q = reduced_form_model.state.Q;
H = reduced_form_model.measurement.H;

smpl = size(Y,2);                             % Sample size.
mm   = size(dr.ghx,2);                        % Number of state variables.
pp   = size(Y,1);                             % Maximum number of observed variables.
qq   = length(Q);
lik  = NaN(smpl,1);
LIK  = 0 ;

if nargin==4
  number_of_particles = 10 ;
  method = 'mc' ;
end
if nargin==5 
  method = 'mc' ;
  if isempty(number_of_particles)
    number_of_particles = 5000 ;
  end    
end
if nargin==6 
  method = 'quadra' ;
  if isempty(grid_size)
    grid_size = 4 ;
  end
end
if isempty(start)
   start = 1; 
end

k2 = dr.kstate(find(dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*M_.endo_nbr;

moy_s_t_1_t_1 = dr.ys ;
tot = M_.endo_nbr ;
Pss_t_1_t_1 = lyapunov_symm(dr.ghx(k2,:),dr.ghu(k2,:)*Q*dr.ghu(k2,:)',1e-12,1e-12);

[V,D] = eig(Pss_t_1_t_1) ;

sqrtP = V*sqrt(D)*V';

dim_var = qq + mm ;
switch method
    case 'mc' 
      %seed  = [ 362436069 ; 521288629 ];
      %randn('state',seed);
      nodes = randn(number_of_particles,dim_var) ;
      weights = 1/number_of_particles ;
      nb_grid = number_of_particles ;
    case 'quadra'
      %{nodes,weights} = integration_smolyak(dim_var,grid_size) ;
      nb_grid = size(nodes,1) ;
    otherwise
      error('undefined method') ; 
end

var_idx = 1:mm;
inn_idx = mm+1:dim_var;
chol_Q = chol(Q);
iorder = 2;
s_t_1_t_1 = zeros(nb_grid,tot) ;


iorder=2;

for t=1:smpl 
  s_t_1_t_1(:,k2) = (repmat(moy_s_t_1_t_1(k2),1,nb_grid) + sqrtP*nodes(:,var_idx)')' ;     
  e_t = (chol_Q*nodes(:,inn_idx)')' ;
  s_t_t_1 = 0 ;
  y_t_t_1 = 0 ;
  Pss_t_t_1 = 0 ;
  Pyy_t_t_1 = 0 ;
  Psy_t_t_1 = 0 ;
  for i=1:nb_grid 
    tout = simult_(s_t_1_t_1(i,:)',dr,e_t(i,:),iorder) ; 
    s_t_t_1 = s_t_t_1 + tout(k2,2) ;
    y_t_t_1 = y_t_t_1 + tout(mf,2) ;
    Pyy_t_t_1  = Pyy_t_t_1  + tout(mf,2)*tout(mf,2)' ; 
    Psy_t_t_1 = Psy_t_t_1 + tout(k2,2)*tout(mf,2)' ;
    Pss_t_t_1 = Pss_t_t_1 + tout(k2,2)*tout(k2,2)' ;
  end
  s_t_t_1 = s_t_t_1/nb_grid ; 
  y_t_t_1 = y_t_t_1/nb_grid ;
  Pyy_t_t_1  = Pyy_t_t_1 + H - y_t_t_1*(y_t_t_1') ;
  Psy_t_t_1  = Psy_t_t_1 - s_t_t_1*(y_t_t_1'); 
  Pss_t_t_1  = Pss_t_t_1 - s_t_t_1*(s_t_t_1'); 
  inv_Pyy_t_t_1 = inv(Pyy_t_t_1) ;
  eta_t = Y(:,t) - y_t_t_1 ;
  Kt = Psy_t_t_1*inv_Pyy_t_t_1 ;
  moy_s_t_1_t_1(k2) = s_t_t_1 + Kt*eta_t ;
  Pss_t_t_1 = Pss_t_t_1 - Kt*Pyy_t_t_1 *(Kt') ;
  [V,D] = eig(Pss_t_1_t_1);
  sqrtP = V*sqrt(D)*V';
  lik(t) = -0.5*log(2*pi)*pp - 0.5*log(det(Pyy_t_t_1)) - 0.5*eta_t'*inv_Pyy_t_t_1*eta_t ;
  if abs(imag(lik(t)))<1e-12  
    lik(t) = real(lik(t));  
  end
end

LIK = -sum(lik(start:end));
