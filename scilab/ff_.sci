function [z]=ff_(y)
z=[];
// 
z = zeros(2,1);

z(1) = y(2)+y(3)-aa*ex_(it_-1,1)*(y(1)^alph)-(1-delt)*y(1)
%v = y(2)^(-gam)-((1+bet)^(-1))*(aa*alph*ex_(it_,1)*(y(3)^(alph-1))+1-delt)*(y(4)^(-gam))
z(2,1) = %v(:);
