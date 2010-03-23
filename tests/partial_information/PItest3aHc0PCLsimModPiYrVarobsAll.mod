% price markup shock (eps_m) AND inflation target (pitarg) modelled as AR1 
% estimating varrho
% cy = 0.614479/0.769365 - obtained from the computed steady state from the nonlinear counterpart
% use usdata1.mat
% this version - 13/07/09


var pi mc mun muc c y n r g a pitarg PI Y RR;
varexo eps_g eps_a eps_e eps_m eps_targ; 

parameters beta xi hc wd sigma gamma rho_g rho_a rho_r rho_targ thetap cy varrho; 
beta = 0.99;
xi = 0.5034;//0.5575;//0.64;
hc = 0.0;//0.7984;//0.86;
wd = 0.40;
gamma = 0.5868;//0.5792;//0.00;
sigma = 4.0897;//2.0034;//3.47 
cy = 0.614479/0.769365; // cy = Cbar/Ybar
rho_g = 0.8325;//0.9132;//0.85;
rho_a = 0.9827;//0.9787;//0.87;
rho_r = 0.3529;//0.4945;//0.76;
thetap = 2.2161;//2.2062;//2.30;
varrho = 0.3853;//0.4432;//0.66;
rho_targ = 0.6133;//0.5807;//0.9;

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
PI = pi;
Y = y;
RR = r;
end;

shocks;
var eps_g; stderr 3.8505;//4.1018;//0.325;
var eps_a; stderr 0.7573;//0.9725;//0.598;
var eps_e; stderr 0.2409;//0.2466;//0.150;
var eps_m; stderr 0.8329;//0.6487;//0.2;
var eps_targ; stderr 0.3978;//0.4577;//0.03;
end;

//steady;
//check;
varobs pi mc mun muc c y n r g a pitarg;
stoch_simul(linear,partial_information,periods=1000,irf=30)pi y r;//pi n c y r;

