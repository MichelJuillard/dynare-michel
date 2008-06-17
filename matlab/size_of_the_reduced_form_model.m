function mega = size_of_the_reduced_form_model(dr)
    names = fieldnames(dr);
    number_of_scalars = 0;
    for field=1:size(names,1)
        number_of_scalars = number_of_scalars + prod(size(getfield(dr,names{field})));
    end
    mega = 8 * number_of_scalars / 1048576 ;