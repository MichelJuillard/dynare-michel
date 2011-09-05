function sp = setSize(sp,n)
%@info:
%! @deftypefn {Function File} {@var{sp} =} setSize (@var{sp}, @var{n})
%! @anchor{@dynTime/setSize}
%! @sp 1
%! Set the size of a dynTime object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! dynTime object instantiated by @ref{dynTime}
%! @item n
%! Positive scalar integer, the number of periods.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! Updated @ref{dynTime} object.
%! @end table
%! @sp 2
%! @strong{Example}
%! @sp 1
%! Let @var{sp} be an object instantiated by @ref{dynTime}, both following syntaxes are equivalent:
%! @sp 1
%! @example
%! sp = setSize(sp,167);
%! @end example
%! or
%! @example
%! sp = sp.setSize(167);
%! @end example
%! @sp 1
%! Note that the second syntax is probably slower than the first one, and should not be used in a large loop.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:
    sp.time = NaN(n,2);
