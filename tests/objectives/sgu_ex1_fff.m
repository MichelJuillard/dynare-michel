
function z=sgu_ex1_fff(y)
z=zeros(3,1);
global ex_ it_ recur_

global alpha beta delta gamma rho 
z(1) = exp(y(2))+exp(y(3)) -((1-delta)*exp(y(3))+exp(y(1))*exp(y(3))^alpha) ...
;
z(2) = exp(y(2))^(-gamma) -(beta*exp(y(2))^(-gamma)*(exp(y(1))*alpha*exp( ...
y(3))^(alpha-1)+1-delta));
z(3) = y(1) -(rho*y(1)+ex_(it_-1,1));
