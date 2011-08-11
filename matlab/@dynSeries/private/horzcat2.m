function a = horzcat2(b,c)

%@info:
%! @deftypefn {Function file} {@var{a} =} horzcat2 (@var{b},@var{c}, ...)
%! @anchor{private/horzcat2}
%! @sp 1
%! Private method of the dynSeries class.
%! @sp 1
%! Merge two Dynare time series objects. This method overloads the horizontal concatenation operator, so that
%! two time series objects can be merged using the following syntax
%!
%!     a = [b, c];
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @item c
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}
%!
%! @strong{This function calls:}
%! @ref{dynSeries}
%!
%! @strong{Remark 1.} It is assumed that the two time series objects have the same frequencies. The two time series objects can cover
%! different time ranges.
%!
%! @end deftypefn
%@eod:

if ~(isa(b,'dynSeries') && isa(c,'dynSeries'))
    error('dynSeries::horzcat: All input arguments have to be Dynare time series objects!')
end

if b.freq ~= c.freq
    error('dynSeries::horzcat: All time series objects must have common frequency!')
else
    a = dynSeries();
    a.freq = b.freq;
end

d_nobs_flag = 0;
if b.nobs ~= c.nobs
    d_nobs_flag = 1;
else
    a.nobs = b.nobs;
end

d_init_flag = 0;
if isequal(b.init,c.init)
    a.init = b.init;
else
    % set a.init equal to min(b.init,c.init)
    if b.init(1)<c.init(1)
        d_init_flag = 1;
        a.init = b.init;
    elseif b.init(1)==c.init(1)
        if b.init(2)<c.init(2)
            d_init_flag = 1;
            a.init = b.init;
        else
            d_init_flag = 2;
            a.init = c.init;
        end
    else
        d_init_flag = 2;
        a.init = c.init;
    end
end

d_last_flag = 0;
if isequal(b.last,c.last)
    a.last = b.last;
else
    % set a.last equal to max(b.last,c.last)
    if b.last(1)<c.last(1)
        d_last_flag = 2;
        a.last = c.last;
    elseif b.last(1)==c.last(1)
        if b.last(2)<c.last(2)
            d_last_flag = 2;
            a.last = c.last;
        else
            d_last_flag = 1;
            a.last = b.last;
        end
    else
        d_last_flag = 1;
        a.last = b.last;
    end
end

a.vobs = b.vobs+c.vobs;
a.name = char(b.name,c.name);
a.tex  = char(b.tex,c.tex);

if ~( d_nobs_flag(1) || d_init_flag(1) || d_last_flag(1) )
    a.time = b.time;
    a.data = [b.data,c.data];
else
    [junk,ib] = setdiff(b.time,c.time,'rows');
    [junk,ic] = setdiff(c.time,b.time,'rows');
    [junk,jb,jc] = intersect(b.time,c.time,'rows');
    a.time = [b.time(ib,:); b.time(jb,:); c.time(ic,:)];
    a.time = sortrows(a.time,[1 2]);
    a.nobs = rows(a.time);
    a.data = NaN(a.nobs,a.vobs);
    [junk,ia,ib] = intersect(a.time,b.time,'rows');
    a.data(ia,1:b.vobs) = b.data;
    [junk,ia,ic] = intersect(a.time,c.time,'rows');
    a.data(ia,b.vobs+(1:c.vobs)) = c.data;
end