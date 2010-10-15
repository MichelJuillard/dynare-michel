// See fs2000.mod in the examples/ directory for details on the model
// This file is the same than fs2000a.mod, except that it is written in non-stationary form
// Notational changes: "m" and "dA" in fs2000a.mod are here called "gM" and "gA"

var gM ${g^M}$ gA ${g^A}$;
trend_var(growth_factor=gA) A;
trend_var(growth_factor=gM) M;
var(deflator=A) k c y;
var(deflator=M(-1)/A) P;
var(deflator=M(-1)) W l d;
var R n;
varexo e_a ${e^A}$ e_m ${e^M}$;

parameters alp $\alpha$ bet $\beta$ gam $\gamma$ mst rho $\rho$ psi $\psi$ del $\delta$;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
gA = exp(gam+e_a);
log(gM) = (1-rho)*log(mst) + rho*log(gM(-1))+e_m;
c+k = k(-1)^alp*(A*n)^(1-alp)+(1-del)*k(-1);
P*c = M;
P/(c(+1)*P(+1))=bet*P(+1)*(alp*k^(alp-1)*(A(+1)*n(+1))^(1-alp)+(1-del))/(c(+2)*P(+2));
(psi/(1-psi))*(c*P/(1-n))=W;
R = P*(1-alp)*k(-1)^alp*A^(1-alp)*n^(-alp)/W;
W = l/n;
M-M(-1)+d = l;
1/(c*P)=bet*R/(c(+1)*P(+1));
y = k(-1)^alp*(A*n)^(1-alp);
end;

initval;
k = 6;
gM = mst;
P = 2.25;
c = 0.45;
W = 4;
R = 1.02;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gA = exp(gam);
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

write_latex_dynamic_model;

stoch_simul(nograph);