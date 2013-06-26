function display(A)
%@info:
%! @deftypefn {Function File} display (@var{A})
%! @anchor{@dynSeries/display}
%! @sp 1
%! Overloads the disp method for the Dynare time series class (@ref{dynSeries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! None
%! @end deftypefn
%@eod:

separator = repmat(' | ',A.nobs+1,1);
vspace = ' ';
TABLE = ' ';
for t=1:A.nobs
    TABLE = char(TABLE, format(A.time(t)));
end
for i = 1:A.vobs
    TABLE = horzcat(TABLE,separator);
    tmp = A.name{i};
    for t=1:A.nobs
        tmp = char(tmp,num2str(A.data(t,i)));
    end
    TABLE = horzcat(TABLE, tmp);
end
disp(vspace)
disp([inputname(1) ' is a dynSeries object:'])
disp(vspace);
disp(TABLE);
disp(vspace);