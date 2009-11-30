// checks for order = 2 that use_k_order = 0 (fs2000k2a) and use_k_order = 1 (this file)
// give the same results

var m P c e W R k d n l gy_obs gp_obs y dA ;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model (use_dll);
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

initval;
m = mst;
P = 2.25;
c = 0.45;
e = 1;
W = 4;
R = 1.02;
k = 6;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gy_obs = exp(gam);
gp_obs = exp(-gam); 
dA = exp(gam);
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

stoch_simul(order=2,use_k_order,periods=1000);

oo1 = load('fs2000k2a_results','oo_');

dr0 = oo1.oo_.dr;
dr = oo_.dr;

if max(max(abs(dr0.ghx - dr.ghx))) > 1e-12;
   disp('error in ghx');
end;
if max(max(abs(dr0.ghu - dr.ghu))) > 1e-12;
   disp('error in ghu');
end;
if max(max(abs(dr0.ghxx - dr.ghxx))) > 1e-12;
   disp('error in ghxx');
end;
if max(max(abs(dr0.ghuu - dr.ghuu))) > 1e-12;
   disp('error in ghuu');
end;
if max(max(abs(dr0.ghxu - dr.ghxu))) > 1e-12;
   disp('error in ghxu');
end;
if max(max(abs(dr0.ghs2 - dr.ghs2))) > 1e-12;
   disp('error in ghs2');
end;

