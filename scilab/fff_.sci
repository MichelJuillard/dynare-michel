function [z]=fff_(y)
z=[];
// 
z = zeros(2,1);
global('ex_','it_','recur_');
 
global('aa','alph','bet','delt','gam');
z(1) = y(1)+y(2)-aa*ex_(it_-1,1)*(y(2)^alph)-(1-delt)*y(2)
%v = y(1)^(-gam)-((1+bet)^(-1))*(aa*alph*ex_(it_,1)*(y(2)^(alph-1))+1-delt)*(y(1)^(-gam))
z(2,1) = %v(:);
