function i = getVarIndex(ts,name)
    switch size(name,1)
      case 0
        error('dynSeries::getVarIndex: Second input argument is empty!');
      case 1
        i = strmatch(deblank(name),ts.name,'exact');
        if isempty(i)
            i = 0;
        end
      otherwise
        i = NaN(size(name,1))
        for j = 1:size(name,1)
            i(j) = strmatch(deblank(name(j,:)),ts.name,'exact');
            if isempty(i)
                i(j) = 0;
            end
        end
    end