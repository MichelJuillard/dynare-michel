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

for j=1:size(x,2);
meany(j)=mean(x(find(~isnan(x(:,j))),j));
stdy(j)=std(x(find(~isnan(x(:,j))),j));
    y(:,j)=(x(:,j)-meany(j))./stdy(j);
end
% end of m-file