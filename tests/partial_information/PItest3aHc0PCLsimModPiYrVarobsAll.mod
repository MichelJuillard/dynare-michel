% inflation target (pitarg) modelled as AR1 
% cy = 0.614479/0.769365 - obtained from the computed steady state from the nonlinear counterpart


var pi mc mun muc c y n r g a pitarg ;
varexo eps_g eps_a eps_e eps_m eps_targ; 

parameters beta xi hc wd sigma gamma rho_g rho_a rho_r rho_targ thetap cy varrho; 
beta = 0.99;
xi = 0.5034;
hc = 0.0;
wd = 0.40;
gamma = 0.5868;
sigma = 4.0897;
cy = 0.614479/0.769365; 
rho_g = 0.8325;
rho_a = 0.9827;
rho_r = 0.3529;
thetap = 2.2161;
varrho = 0.3853;
rho_targ = 0.6133;

model(linear);
pi = (beta/(1+beta*gamma))*pi(+1)+(gamma/(1+beta*gamma))*pi(-1)+(((1-beta*xi)*(1-xi))/((1+beta*gamma)*xi))*(mc+eps_m);
mc = mun-muc-a;
mun = (c-hc*c(-1))/(1-hc)+wd*n/(1-wd)+muc;
muc = ((1-varrho)*(1-sigma)-1)*(c-hc*c(-1))/(1-hc)-wd*varrho*(1-sigma)*n/(1-wd);
muc(+1) = muc-(r-pi(+1));
y = cy*c+(1-cy)*g;
n = y-a;
r = rho_r*r(-1)+thetap*(1-rho_r)*(pi(+1)-rho_targ*pitarg)+eps_e;
g = rho_g*g(-1)+eps_g;
a = rho_a*a(-1)+eps_a;
pitarg = rho_targ*pitarg(-1)+eps_targ;
end;

shocks;
var eps_g; stderr 3.8505;
var eps_a; stderr 0.7573;
var eps_e; stderr 0.2409;
var eps_m; stderr 0.8329;
var eps_targ; stderr 0.3978;
end;

varobs pi mc mun muc c y n r g a pitarg;
stoch_simul(partial_information,irf=30)pi y r;//pi n c y r;

