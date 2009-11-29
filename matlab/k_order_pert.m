function [dr,info] = k_order_pert(dr,M,options,oo)
    
    info = 0;
    
    M.var_order_endo_names = M.endo_names(dr.order_var,:);

    order = options.order;
    
    switch(order)
      case 1
        g_1 = k_order_perturbation(dr,0,M,options, oo , ['.' ...
                            mexext]);
      case 2
        [g_0, g_1, g_2] = k_order_perturbation(dr,0,M,options, oo , ['.' ...
                            mexext]);
      case 3
        [g_0, g_1, g_2, g_3] = k_order_perturbation(dr,0,M,options, oo , ['.' ...
                            mexext]);
      otherwise
        error('order > 3 isn''t implemented')
    end
    
    npred = dr.npred;
    
    dr.ghx = g_1(:,1:npred);
    dr.ghu = g_1(:,npred+1:end);
    
    if order > 1
        dr.ghs2 = 2*g_0;
        endo_nbr = M.endo_nbr;
        exo_nbr = M.exo_nbr;
        s0 = 0;
        s1 = 0;
        ghxx=zeros(endo_nbr, npred^2);
        ghxu=zeros(endo_nbr, npred*exo_nbr);
        ghuu=zeros(endo_nbr, exo_nbr^2);
        for i=1:size(g_2,2)
            if s0 < npred & s1 < npred
                ghxx(:,s0*npred+s1+1) = 2*g_2(:,i);
                if s1 > s0
                    ghxx(:,s1*npred+s0+1) = 2*g_2(:,i);
                end
            elseif s0 < npred & s1 < npred+exo_nbr 
                ghxu(:,(s0*exo_nbr+s1-npred+1)) = 2*g_2(:,i);
            elseif s0 < npred+exo_nbr & s1 < npred+exo_nbr
                ghuu(:,(s0-npred)*exo_nbr+s1-npred +1) = 2*g_2(:,i);
                if s1 > s0
                    ghuu(:,(s1-npred)*exo_nbr+s0-npred+1) = 2*g_2(:,i);
                end
            else
                error('dr1:k_order_perturbation:g_2','Unaccounted columns in g_2');
            end
            s1 = s1+1;
            if s1 == npred+exo_nbr
                s0 = s0+1;
                s1 = s0; 
            end
        end % for loop
        dr.ghxx = ghxx;
        dr.ghxu = ghxu;
        dr.ghuu = ghuu;
    end
    
    if order > 2
        dr.g_0 = g_0;
        dr.g_1 = g_1;
        dr.g_2 = g_2;
        dr.g_3 = g_3;
    end