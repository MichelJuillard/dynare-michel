function h = cumplot(x);
%function h =cumplot(x)
% Copyright (C) 2005 Marco Ratto


n=length(x);
x=[-inf; sort(x); Inf];
y=[0:n n]./n;
h0 = stairs(x,y);
grid on,

if nargout,
    h=h0;
end
