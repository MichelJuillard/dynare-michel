function [dr,info]=AIM_first_order_solver(jacobia,M,dr,qz_criterium,nd)
    
    info = 0;
    
    [dr,aimcode]=dynAIMsolver1(jacobia,M,dr);

    if aimcode ~=1
        info(1) = aimcode;
        info(2) = 1.0e+8;
        return
    end
    [A,B] =transition_matrix(dr);
    dr.eigval = eig(A);
    nd = size(dr.kstate,1);
    nba = nd-sum( abs(dr.eigval) < qz_criterium );

    nyf = dr.nfwrd+dr.nboth;

    if nba ~= nyf
        temp = sort(abs(dr.eigval));
        if nba > nyf
            temp = temp(nd-nba+1:nd-nyf)-1-qz_criterium;
            info(1) = 3;
        elseif nba < nyf;
            temp = temp(nd-nyf+1:nd-nba)-1-qz_criterium;
            info(1) = 4;
        end
        info(2) = temp'*temp;
        return
    end

