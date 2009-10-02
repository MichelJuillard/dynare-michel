function ys1 = add_auxiliary_variables_to_steadystate(ys,aux_vars)
    n = length(aux_vars);
    ys1 = [ys;zeros(n,1)];
    k = size(ys,1)+1;
    for i=1:n
        ys1(k) = ys(aux_vars(i).orig_endo_index);
        k = k+1;
    end
        