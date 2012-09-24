function c = step_length_correction(x,scale,i)
    if isempty(scale)
        c = 10^round(log10(abs(x)));
    else
        c = scale(i);
    end