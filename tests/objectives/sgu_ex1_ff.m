
function z=sgu_ex1_ff(y)
z=zeros(3,1);
global ex_ it_ recur_

global alpha beta delta gamma rho 
z(1) = exp(y(4))+exp(y(5)) -((1-delta)*exp(y(2))+exp(y(3))*exp(y(2))^alpha) ...
;
z(2) = exp(y(4))^(-gamma) -(beta*exp(y(7))^(-gamma)*(exp(y(6))*alpha*exp( ...
y(5))^(alpha-1)+1-delta));
z(3) = y(3) -(rho*y(1)+ex_(it_-1,1));
