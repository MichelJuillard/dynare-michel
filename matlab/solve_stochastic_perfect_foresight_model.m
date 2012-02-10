function [flag,endo_simul,err] = solve_perfect_foresight_model(endo_simul,exo_simul,pfm)

    flag = 0;
    err = 0;
    stop = 0;

    model_dynamic = pfm.dynamic_model;

    Y = endo_simul(:);

    if pfm.verbose
        disp (['-----------------------------------------------------']) ;
        disp (['MODEL SIMULATION :']) ;
        fprintf('\n') ;
    end

    z = Y(find(pfm.lead_lag_incidence'));
    [d1,jacobian] = model_dynamic(z,exo_simul,pfm.params,pfm.steady_state,2);

    A = sparse([],[],[],pfm.periods*pfm.ny,pfm.periods*pfm.ny,pfm.periods*nnz(jacobian));
    res = zeros(pfm.periods*pfm.ny,1);

    h1 = clock;
    for iter = 1:pfm.maxit_
        h2 = clock;
        i_rows = 1:pfm.ny;
        i_cols = find(pfm.lead_lag_incidence');
        i_cols_A = i_cols;
        for it = 2:(pfm.periods+1)
            [d1,jacobian] = model_dynamic(Y(i_cols),exo_simul,pfm.params,pfm.steady_state,it);
            if it == 2
                A(i_rows,pfm.i_cols_A1) = jacobian(:,pfm.i_cols_1);
            elseif it == pfm.periods+1
                A(i_rows,i_cols_A(pfm.i_cols_T)) = jacobian(:,pfm.i_cols_T);
            else
                A(i_rows,i_cols_A) = jacobian(:,pfm.i_cols_j);
            end
            res(i_rows) = d1;
            i_rows = i_rows + pfm.ny;
            i_cols = i_cols + pfm.ny;
            if it > 2
                i_cols_A = i_cols_A + pfm.ny;
            end
        end
        err = max(abs(res));
        if err < pfm.tolerance
            stop = 1 ;
            if pfm.verbose
                fprintf('\n') ;
                disp([' Total time of simulation        :' num2str(etime(clock,h1))]) ;
                fprintf('\n') ;
                disp([' Convergency obtained.']) ;
                fprintf('\n') ;
            end
            flag = 0;% Convergency obtained.
            endo_simul = reshape(Y,pfm.ny,pfm.periods+2);
            break
        end
        dy = -A\res;
        Y(pfm.i_upd) =   Y(pfm.i_upd) + dy;
    end

    if ~stop
        if pfm.verbose
            fprintf('\n') ;
            disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
            fprintf('\n') ;
            disp(['WARNING : maximum number of iterations is reached (modify options_.maxit_).']) ;
            fprintf('\n') ;
        end
        flag = 1;% more iterations are needed.
        endo_simul = 1;
    end
    if pfm.verbose
        disp (['-----------------------------------------------------']) ;
    end