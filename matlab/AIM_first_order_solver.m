function [dr,info]=AIM_first_order_solver(jacobia,M_,dr)
    try
        [dr,aimcode]=dynAIMsolver1(jacobia_,M_,dr);

        % reuse some of the bypassed code and tests that may be needed 
        
        if aimcode ~=1
            info(1) = aimcode;
            info(2) = 1.0e+8;
            return
        end
        [A,B] =transition_matrix(dr);
        dr.eigval = eig(A);
        nba = nd-sdim;

        nyf = sum(kstate(:,2) > M_.maximum_endo_lag+1);

        if nba ~= nyf
            temp = sort(abs(dr.eigval));
            if nba > nyf
                temp = temp(nd-nba+1:nd-nyf)-1-options_.qz_criterium;
                info(1) = 3;
            elseif nba < nyf;
                temp = temp(nd-nyf+1:nd-nba)-1-options_.qz_criterium;
                info(1) = 4;
            end
            info(2) = temp'*temp;
            return
        end
    catch
        disp(lasterror.message)
        error('Problem with AIM solver - Try to remove the "aim_solver" option')
    end
