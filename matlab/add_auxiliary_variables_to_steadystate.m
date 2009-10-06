function ys1 = add_auxiliary_variables_to_steadystate(ys,aux_vars,fname, ...
                                                      exo_steady_state, exo_det_steady_state,params)
    n = length(aux_vars);
    ys1 = [ys;zeros(n,1)];
    k = size(ys,1)+1;
    aux_lead_nbr = 0;
    for i=1:n
        if aux_vars(i).type == 1
            ys1(k) = ys(aux_vars(i).orig_endo_index);
        elseif aux_vars(i).type == 0
            aux_lead_nbr = aux_lead_nbr + 1;
        end
        k = k+1;
    end
    
    for i=1:aux_lead_nbr + 1;
        res = feval([fname '_static'],ys1,...
			 [exo_steady_state; ...
                      exo_det_steady_state],params);
        for j=1:n
            if aux_vars(j).type == 0
                el = aux_vars(j).endo_index;
                ys1(el) == res(el);
            end
        end
    end
    