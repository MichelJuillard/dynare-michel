function sim1()
% function sim1
% Performs deterministic simulations with lead or lag on one period.
% Uses sparse matrices.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 1996-2012 Dynare Team
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

global M_ options_ oo_

lead_lag_incidence = M_.lead_lag_incidence;

ny = M_.endo_nbr;

max_lag = M_.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
iyp = find(lead_lag_incidence(1,:)>0) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
iy0 = find(lead_lag_incidence(2,:)>0) ;
nyf = nnz(lead_lag_incidence(3,:)) ;
iyf = find(lead_lag_incidence(3,:)>0) ;

nd = nyp+ny0+nyf;
nrc = nyf+1 ;
isp = [1:nyp] ;
is = [nyp+1:ny+nyp] ;
isf = iyf+nyp ;
isf1 = [nyp+ny+1:nyf+nyp+ny+1] ;
stop = 0 ;
iz = [1:ny+nyp+nyf];

periods = options_.periods
steady_state = oo_.steady_state;
params = M_.params;
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;
i_cols_1 = nonzeros(lead_lag_incidence(2:3,:)');
i_cols_A1 = find(lead_lag_incidence(2:3,:)');
i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
i_cols_j = 1:nd;
i_upd = ny+(1:periods*ny);

Y = endo_simul(:);

disp (['-----------------------------------------------------']) ;
disp (['MODEL SIMULATION :']) ;
fprintf('\n') ;


model_dynamic = str2func([M_.fname,'_dynamic']);
z = Y(find(lead_lag_incidence'));
[d1,jacobian] = model_dynamic(z,oo_.exo_simul, params, ...
                              steady_state,2);

A = sparse([],[],[],periods*ny,periods*ny,periods*nnz(jacobian));
res = zeros(periods*ny,1);

    
h1 = clock ;
for iter = 1:options_.maxit_
    h2 = clock ;
    
    i_rows = 1:ny;
    i_cols = find(lead_lag_incidence');
    i_cols_A = i_cols;
    
    for it = 2:(periods+1)

        [d1,jacobian] = model_dynamic(Y(i_cols),exo_simul, params, ...
                                      steady_state,it);
        if it == 2
            A(i_rows,i_cols_A1) = jacobian(:,i_cols_1);
        elseif it == periods+1
            A(i_rows,i_cols_A(i_cols_T)) = jacobian(:,i_cols_T);
        else
            A(i_rows,i_cols_A) = jacobian(:,i_cols_j);
        end

        res(i_rows) = d1;
        
        i_rows = i_rows + ny;
        i_cols = i_cols + ny;
        if it > 2
            i_cols_A = i_cols_A + ny;
        end
    end
        
    err = max(abs(res));
    
    if err < options_.dynatol.f
        stop = 1 ;
        fprintf('\n') ;
        disp([' Total time of simulation        :' num2str(etime(clock,h1))]) ;
        fprintf('\n') ;
        disp([' Convergency obtained.']) ;
        fprintf('\n') ;
        oo_.deterministic_simulation.status = 1;% Convergency obtained.
        oo_.deterministic_simulation.error = err;
        oo_.deterministic_simulation.iterations = iter;
        oo_.endo_simul = reshape(Y,ny,periods+2);
        break
    end

    dy = -A\res;
    
    Y(i_upd) =   Y(i_upd) + dy;

end


if ~stop
    fprintf('\n') ;
    disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
    fprintf('\n') ;
    disp(['WARNING : maximum number of iterations is reached (modify options_.maxit_).']) ;
    fprintf('\n') ;
    oo_.deterministic_simulation.status = 0;% more iterations are needed.
    oo_.deterministic_simulation.error = err;
    oo_.deterministic_simulation.errors = c/abs(err);    
    oo_.deterministic_simulation.iterations = options_.maxit_;
end
disp (['-----------------------------------------------------']) ;

