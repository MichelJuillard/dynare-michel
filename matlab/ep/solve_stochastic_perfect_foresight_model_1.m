function [flag,endo_simul,err] = solve_stochastic_perfect_foresight_model_1(endo_simul,exo_simul,pfm,nnodes,order)

% Copyright (C) 2012-2013 Dynare Team
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

    flag = 0;
    err = 0;
    stop = 0;

    params = pfm.params;
    steady_state = pfm.steady_state;
    ny = pfm.ny;
    periods = pfm.periods;
    dynamic_model = pfm.dynamic_model;
    lead_lag_incidence = pfm.lead_lag_incidence;
    nyp = pfm.nyp;
    nyf = pfm.nyf;
    i_cols_1 = pfm.i_cols_1;
    i_cols_A1 = pfm.i_cols_A1;
    i_cols_j = pfm.i_cols_j;
    i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
    hybrid_order = pfm.hybrid_order;
    dr = pfm.dr;
    
    maxit = pfm.maxit_;
    tolerance = pfm.tolerance;
    verbose = pfm.verbose;

    number_of_shocks = size(exo_simul,2);

    [nodes,weights] = gauss_hermite_weights_and_nodes(nnodes);
    
    % make sure that there is a node equal to zero
    % and permute nodes and weights to have zero first
    k = find(abs(nodes) < 1e-12);
    if ~isempty(k)
        nodes = [nodes(k); nodes(1:k-1); nodes(k+1:end)];
        weights = [weights(k); weights(1:k-1); weights(k+1:end)];
    else
        error('there is no nodes equal to zero')
    end
   
    if number_of_shocks>1
        nodes = repmat(nodes,1,number_of_shocks)*chol(pfm.Sigma);
        % to be fixed for Sigma ~= I
        for i=1:number_of_shocks
            rr(i) = {nodes(:,i)};
            ww(i) = {weights};
        end
        nodes = cartesian_product_of_sets(rr{:});
        weights = prod(cartesian_product_of_sets(ww{:}),2);
        nnodes = nnodes^number_of_shocks;
    else
        nodes = nodes*sqrt(pfm.Sigma);
    end
    
    if hybrid_order > 0
        if hybrid_order == 2
            h_correction = 0.5*dr.ghs2(dr.inv_order_var);
        end
    end

    if verbose
        disp ([' -----------------------------------------------------']);
        disp (['MODEL SIMULATION :']);
        fprintf('\n');
    end

    z = endo_simul(find(lead_lag_incidence'));
    [d1,jacobian] = dynamic_model(z,exo_simul,params,steady_state,2);

    % Each column of Y represents a different world
    % The upper right cells are unused
    % The first row block is ny x 1
    % The second row block is ny x nnodes
    % The third row block is ny x nnodes^2
    % and so on until size ny x nnodes^order
    world_nbr = 1+(nnodes-1)*order;
    Y = repmat(endo_simul(:),1,world_nbr);

    % The columns of A map the elements of Y such that
    % each block of Y with ny rows are unfolded column wise
    dimension = ny*(order+(nnodes-1)*(order-1)*order/2+(periods-order)*world_nbr);
    if order == 0
        i_upd_r = (1:ny*periods);
        i_upd_y = i_upd_r + ny;
    else
        i_upd_r = zeros(dimension,1);
        i_upd_y = i_upd_r;
        i_upd_r(1:ny) = (1:ny);
        i_upd_y(1:ny) = ny+(1:ny);
        i1 = ny+1;
        i2 = 2*ny;
        n1 = ny+1;
        n2 = 2*ny;
        for i=2:periods
            k = n1:n2;
            for j=1:(1+(nnodes-1)*min(i-1,order))
                i_upd_r(i1:i2) = (n1:n2)+(j-1)*ny*periods;
                i_upd_y(i1:i2) = (n1:n2)+ny+(j-1)*ny*(periods+2);
                i1 = i2+1;
                i2 = i2+ny;
            end
            n1 = n2+1;
            n2 = n2+ny;
        end
    end
    icA = [find(lead_lag_incidence(1,:)) find(lead_lag_incidence(2,:))+world_nbr*ny ...
           find(lead_lag_incidence(3,:))+2*world_nbr*ny]';
    h1 = clock;
    for iter = 1:maxit
        h2 = clock;
        A1 = sparse([],[],[],ny*(order+(nnodes-1)*(order-1)*order/2),dimension,(order+1)*world_nbr*nnz(jacobian));
        res = zeros(ny,periods,world_nbr);
        i_rows = 1:ny;
        i_cols = find(lead_lag_incidence');
        i_cols_p = i_cols(1:nyp);
        i_cols_s = i_cols(nyp+(1:ny));
        i_cols_f = i_cols(nyp+ny+(1:nyf));
        i_cols_A = i_cols;
        i_cols_Ap0 = i_cols_p;
        i_cols_As = i_cols_s;
        i_cols_Af0 = i_cols_f - ny;
        i_hc = i_cols_f - 2*ny;
        for i = 1:order+1
            i_w_p = 1;
            for j = 1:(1+(nnodes-1)*(i-1))
                innovation = exo_simul;
                if i <= order && j == 1
                    % first world, integrating future shocks
                    for k=1:nnodes
                        if i == 2
                            i_cols_Ap = i_cols_Ap0;
                        elseif i > 2
                            i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes- ...
                                                              1)*(i-2)*(i-3)/2);
                        end
                        if k == 1
                            i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                         2);
                        else
                            i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                         2+(i-1)*(nnodes-1) ...
                                                         + k - 1);
                        end
                        if i > 1
                            innovation(i+1,:) = nodes(k,:);
                        end
                        if k == 1
                            k1 = 1;
                        else
                            k1 = (nnodes-1)*(i-1)+k;
                        end
                        if hybrid_order == 2 && (k > 1 || i == order) 
                            y = [Y(i_cols_p,1);
                                 Y(i_cols_s,1);
                                 Y(i_cols_f,k1)+h_correction(i_hc)];
                        else
                            y = [Y(i_cols_p,1);
                                 Y(i_cols_s,1);
                                 Y(i_cols_f,k1)];
                        end
                        [d1,jacobian] = dynamic_model(y,innovation,params,steady_state,i+1);
                        if i == 1
                            % in first period we don't keep track of
                            % predetermined variables
                            i_cols_A = [i_cols_As - ny; i_cols_Af];
                            A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(:,i_cols_1);
                        else
                            i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                            A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(:,i_cols_j);
                        end
                        res(:,i,1) = res(:,i,1)+weights(k)*d1;
                    end
                elseif j > 1 + (nnodes-1)*(i-2)
                    % new world, using previous state of world 1 and hit
                    % by shocks from nodes
                    i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2);
                    i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
                    k = j - (nnodes-1)*(i-2);
                    innovation(i+1,:) = nodes(k,:);
                    y = [Y(i_cols_p,1);
                         Y(i_cols_s,j);
                         Y(i_cols_f,j)];
                    [d1,jacobian] = dynamic_model(y,innovation,params,steady_state,i+1);
                    i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                    A1(i_rows,i_cols_A) = jacobian(:,i_cols_j);
                    res(:,i,j) = d1;
                    i_cols_Af = i_cols_Af + ny;
                else
                    % existing worlds other than 1
                    i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2+j-1);
                    i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
                    y = [Y(i_cols_p,j);
                         Y(i_cols_s,j);
                         Y(i_cols_f,j)];
                    [d1,jacobian] = dynamic_model(y,innovation,params,steady_state,i+1);
                    i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                    A1(i_rows,i_cols_A) = jacobian(:,i_cols_j);
                    res(:,i,j) = d1;
                    i_cols_Af = i_cols_Af + ny;
                end
                i_rows = i_rows + ny;
                if mod(j,nnodes) == 0
                    i_w_p = i_w_p + 1;
                end
                if i > 1
                    i_cols_As = i_cols_As + ny;
                end
            end
            i_cols_p = i_cols_p + ny;
            i_cols_s = i_cols_s + ny;
            i_cols_f = i_cols_f + ny;
        end
        nzA = cell(periods,world_nbr);
        parfor j=1:world_nbr
            i_rows_y = find(lead_lag_incidence')+(order+1)*ny;
            offset_c = ny*(order+(nnodes-1)*(order-1)*order/2+j-1);
            offset_r = (j-1)*ny;
            for i=order+2:periods
                [d1,jacobian] = dynamic_model(Y(i_rows_y,j), ...
                                              exo_simul,params, ...
                                              steady_state,i+1);
                if i == periods
                    [ir,ic,v] = find(jacobian(:,i_cols_T));
                else
                    [ir,ic,v] = find(jacobian(:,i_cols_j));
                end
                nzA{i,j} = [offset_r+ir,offset_c+icA(ic), v]';
                res(:,i,j) = d1;
                i_rows_y = i_rows_y + ny;
                offset_c = offset_c + world_nbr*ny;
                offset_r = offset_r + world_nbr*ny;
            end
        end
        err = max(abs(res(i_upd_r)));
        if err < tolerance
            stop = 1;
            if verbose
                fprintf('\n') ;
                disp([' Total time of simulation        :' num2str(etime(clock,h1))]) ;
                fprintf('\n') ;
                disp([' Convergency obtained.']) ;
                fprintf('\n') ;
            end
            flag = 0;% Convergency obtained.
            endo_simul = reshape(Y(:,1),ny,periods+2);%Y(ny+(1:ny),1);
                                                      %            figure;plot(Y(16:ny:(periods+2)*ny,:))
                                                      %            pause
            break
        end
        A2 = [nzA{:}]';
        A = [A1; sparse(A2(:,1),A2(:,2),A2(:,3),ny*(periods-order-1)*world_nbr,dimension)];
        dy = -A\res(i_upd_r);
        Y(i_upd_y) =   Y(i_upd_y) + dy;
    end

    if ~stop
        if verbose
            fprintf('\n') ;
            disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
            fprintf('\n') ;
            disp(['WARNING : maximum number of iterations is reached (modify options_.simul.maxit).']) ;
            fprintf('\n') ;
        end
        flag = 1;% more iterations are needed.
        endo_simul = 1;
    end
    if verbose
        disp (['-----------------------------------------------------']) ;
    end