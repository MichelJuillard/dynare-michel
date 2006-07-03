function y=nonzeros(x);
x=x(:);
k = find(x);
y=x(k);
