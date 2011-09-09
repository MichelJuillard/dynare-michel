function sp = setTime(sp,i,date)
    if nargin==3
        sp.time(i,:) = date;
    elseif nargin==2
        if isa(i,'dynTime')
            sp.time=i.time;
        else
            sp.time=i;
        end
    end