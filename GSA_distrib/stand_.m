function [y, meany, stdy] = stand_(x)
% STAND_  Standardise a matrix by columns
%
% [x,my,sy]=stand(y)
%
% y: Time series (column matrix)
%
% x: standardised equivalent of y
% my: Vector of mean values for each column of y
% sy: Vector of standard deviations for each column of y
%
% Copyright (c) 2006 by JRC, European Commission, United Kingdom
% Author : Marco Ratto


if nargin==0,
    return;
end

meany=mean(x);
stdy=std(x);
for j=1:size(x,2);
    y(:,j)=(x(:,j)-meany(j))./stdy(j);
end
% end of m-file