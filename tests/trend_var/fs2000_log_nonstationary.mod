/*
 * This file is a modified version of 'fs2000.mod'.
 *
 * The difference is that, here, the equations are written in non-stationary form,
 * all variables are taken in logs, and Dynare automatically does the detrending.
 *
 * Also note that "m" and "dA" in 'fs2000.mod' are here called "gM" and "gA"
 */

/*
 * Copyright (C) 2004-2013 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

var gM gA;
log_trend_var(log_growth_factor=gA) A;
log_trend_var(log_growth_factor=gM) M;
var(log_deflator=A) k c y;
var(log_deflator=M(-1)-A) P;
var(log_deflator=M(-1)) W l d;
var R n;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
gA = gam+e_a;
gM = (1-rho)*log(mst) + rho*gM(-1)+e_m;
exp(c)+exp(k) = exp(k(-1))^alp*(exp(A)*exp(n))^(1-alp)+(1-del)*exp(k(-1));
P+c = M;
P-(c(+1)+P(+1))=log(bet)+P(+1)+log(alp*exp(k)^(alp-1)*(exp(A(+1)+n(+1)))^(1-alp)+(1-del))-(c(+2)+P(+2));
log(psi/(1-psi))+(c+P-log(1-exp(n)))=W;
R = P+log(1-alp)+alp*k(-1)+(1-alp)*A+(-alp)*n-W;
W = l-n;
exp(M)-exp(M(-1))+exp(d) = exp(l);
-(c+P)=log(bet)+R-(c(+1)+P(+1));
y = alp*k(-1)+(1-alp)*(A+n);
end;

initval;
k = log(6);
gM = log(mst);
P = log(2.25);
c = log(0.45);
W = log(4);
R = log(1.02);
d = log(0.85);
n = log(0.19);
l = log(0.86);
y = log(0.6);
gA = gam;
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

stoch_simul;
