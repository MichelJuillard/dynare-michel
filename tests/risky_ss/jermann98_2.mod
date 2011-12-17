// This is the Ramsey model with adjustment costs.  Jermann(1998),JME 41, pages 257-275
// Olaf Weeken
// Bank of England, 13 June, 2005
// modified January 20, 2006 by Michel Juillard

//---------------------------------------------------------------------
// 1. Variable declaration
//---------------------------------------------------------------------

var c, d, erp1, i, k, m1, r1, rf1, w, y, z, mu; 
varexo ez;                          

//---------------------------------------------------------------------
// 2. Parameter declaration and calibration
//---------------------------------------------------------------------

parameters alf, chihab, xi, delt, tau, g, rho, zbar, a1, a2, betstar, bet;

alf        = 0.36;    // capital share in production function
//chihab     = 0.819;   // habit formation parameter
chihab     = 0.98;   // habit formation parameter
xi         = 1/4.3;   // capital adjustment cost parameter
delt       = 0.025;   // quarterly deprecition rate
g          = 1.005;   //quarterly growth rate (note zero growth =>g=1)
tau        = 5;       // curvature parameter with respect to c
rho        = 0.95;    // AR(1) parameter for technology shock

a1         = (g-1+delt)^(1/xi);             
a2         = (g-1+delt)-(((g-1+delt)^(1/xi))/(1-(1/xi)))*((g-1+delt)^(1-(1/xi))); 
betstar    = g/1.011138;
bet        = betstar/(g^(1-tau));             

//---------------------------------------------------------------------
// 3. Model declaration
//---------------------------------------------------------------------

model;  
g*k  = (1-delt)*k(-1) + ((a1/(1-1/xi))*(g*i/k(-1))^(1-1/xi)+a2)*k(-1);
d    = y - w - i; 
w    = (1-alf)*y;
y    = z*g^(-alf)*k(-1)^alf;
c    = w + d; 
mu   = ((c-chihab*c(-1)/g)^(-tau)-chihab*bet*(c(+1)*g-chihab*c)^(-tau))/1e4;
mu   = (betstar/g)*mu(+1)*(a1*(g*i/k(-1))^(-1/xi))*(alf*z(+1)*g^(1-alf)*
       (k^(alf-1))+((1-delt+(a1/(1-1/xi))*(g*i(+1)/k)^(1-1/xi)+a2))/
       (a1*(g*i(+1)/k)^(-1/xi))-g*i(+1)/k);
log(z) = rho*log(z(-1)) + ez;

m1   = (betstar/g)*mu(+1)/mu;
rf1  = 1/m1;
r1   = (a1*(g*i/k(-1))^(-1/xi))*(alf*z(+1)*g^(1-alf)*(k^(alf-1))+
       (1-delt+(a1/(1-1/xi))*(g*i(+1)/k)^(1-1/xi)+a2)/
       (a1*(g*i(+1)/k)^(-1/xi))-g*i(+1)/k);
erp1 = r1 - rf1;

end;

//---------------------------------------------------------------------
// 4. Initial values and steady state
//---------------------------------------------------------------------

initval;
m1     = betstar/g;
rf1    = (1/m1);
r1     = (1/m1);
erp1   = r1-rf1;

z      = 1;
k      = (((g/betstar)-(1-delt))/(alf*g^(1-alf)))^(1/(alf-1));
y      = (g^(-alf))*k^alf;
w      = (1-alf)*y;
i      = (1-(1/g)*(1-delt))*k;
d      = y - w - i;
c      = w + d;

mu     = (((c-(chihab*c/g))^(-tau))-chihab*bet*((c*g-chihab*c)^(-tau)))/1e4;

ez     = 0;
end;

resid(1);

steady;                      

//---------------------------------------------------------------------
// 5. Shock declaration  
//                       
//---------------------------------------------------------------------

for i=1:100:101;
s = i/10000;
shocks;
var ez; stderr s;  
end;

options_.risky_steadystate = 1;
stoch_simul (order=2,irf=0,noprint) erp1, rf1, m1, r1, y, z, c, d, mu, k;
oo_.steady_state = oo_.dr.ys;
end