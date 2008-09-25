{\rtf1\mac\ansicpg10000\cocoartf824\cocoasubrtf420
{\fonttbl\f0\fswiss\fcharset77 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww9000\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\ql\qnatural\pardirnatural

\f0\fs24 \cf0 // example 1 from Collard's guide to Dynare\
var y, c, k, a, h, b;\
varexo e,u;\
\
parameters beta, rho, beta, alpha, delta, theta, psi, tau;\
\
alpha = 0.36;\
rho   = 0.95;\
tau   = 0.025;\
beta  = 0.99;\
delta = 0.025;\
psi   = 0;\
theta = 2.95;\
\
phi   = 0.1;\
\
model;\
c*theta*h^(1+psi)=(1-alpha)*y;\
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))\
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));\
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));\
k = exp(b)*(y-c)+(1-delta)*k(-1);\
a = rho*a(-1)+tau*b(-1) + e;\
b = tau*a(-1)+rho*b(-1) + u;\
end;\
\
initval;\
y = 1.08068253095672;\
c = 0.80359242014163;\
h = 0.29175631001732;\
k = 5;\
a = 0;\
b = 0;\
e = 0;\
u = 0;\
end;\
\
shocks;\
var e; stderr 0.009;\
var u; stderr 0.009;\
var e, u = phi*0.009*0.009;\
end;\
\
stoch_simul(periods=2100);\
}