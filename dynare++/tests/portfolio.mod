var DOTQ Q1 Q2 X1 X2 C D1 D2;

varexo E_D1 E_D2;

parameters beta, r1, r2, gamma, d, rho1, rho2;

beta = 0.95;
r1 = 0.2;
r2 = 0.05;

gamma = 0.78;
d = 0.10;

rho1 = 0.8;
rho2 = 0.2;

model;
C + X1 + X2 = D1*Q1 + D2*Q2;
Q1+Q2 = 1;
C^(-gamma)/(1-2*r1*X1) = beta*DOTQ(+1)^(-gamma)*C(+1)^(-gamma)/(1-2*r1*X1(+1))*(D1(+1)*(1-2*r1*X1(+1))+1);
C^(-gamma)/(1-2*r2*X2) = beta*DOTQ(+1)^(-gamma)*C(+1)^(-gamma)/(1-2*r2*X2(+1))*(D2(+1)*(1-2*r2*X2(+1))+1);
DOTQ*Q1 = Q1(-1) + X1(-1) - r1*X1(-1)^2;  
DOTQ*Q2 = Q2(-1) + X2(-1) - r2*X2(-1)^2;

D1/d = D1(-1)^rho1/(d^rho1)*exp(E_D1);
D2/d = D2(-1)^rho2/(d^rho2)*exp(E_D2);

/*
D1-d = rho1*(D1(-1)-d) + E_D1;
D2-d = rho2*(D2(-1)-d) + E_D2;
*/
end;

initval;
C    		 =0.0441234;
D1   		 =0.1000000000000;
D2   		 =0.1000000000000;

DOTQ 		 =1.05567;
Q1   		 =0.333333;
Q2   		 =0.666667;

X1   		 =0.0186255;
X2   		 =0.0372511;
end;

vcov = [
0.04 0;
0 0.01
];

order=5;
