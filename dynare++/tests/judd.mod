var C K1 K2 L S1 S2 THETA V V1 V2;

varexo KSI;

parameters theta_ss lambda delta1 delta2 alpha1 alpha2 eta beta gamma depr1 depr2;

theta_ss=1;
lambda=0.8;
delta1=0.1;
delta2=0.05;
alpha1=0.3;
alpha2=0.15;
eta=3;
beta=0.95;
gamma=0.5;
depr1=0.1;
depr2=0.05;

model;
C = THETA*K1^alpha1*K2^alpha2*L^(1-alpha1-alpha2)-S1*K1-S2*K2;
K1 = (1-depr1+(1-0.5*delta1*S1)*S1)*K1(-1);
K2 = (1-depr2+(1-0.5*delta2*S2)*S2)*K2(-1);
THETA = THETA(-1)^lambda*theta_ss^(1-lambda)*exp(KSI);
/*
THETA = THETA(-1)*lambda+theta_ss*(1-lambda)+KSI;
*/
C^(-gamma)*THETA*K1^alpha1*K2^alpha2*L^(-alpha1-alpha2)*(1-alpha1-alpha2) = L^eta;
C^(-gamma) = beta*V1(+1)*(1-delta1*S1); 
C^(-gamma) = beta*V2(+1)*(1-delta2*S2);
V1 = C^(-gamma)*(alpha1*THETA*K1^(alpha1-1)*K2^alpha2*L^(1-alpha1-alpha2)-S1)+beta*V1(+1)*(1-depr1+(1-0.5*delta1*S1)*S1); 
V2 = C^(-gamma)*(alpha2*THETA*K1^alpha1*K2^(alpha2-1)*L^(1-alpha1-alpha2)-S2)+beta*V2(+1)*(1-depr2+(1-0.5*delta2*S2)*S2);
V = (C^(1-gamma)/(1-gamma)-L^(1+eta)/(1+eta)) + beta*V(+1);
end;

initval;
C=   1.33341818203972;
K1=   3.80023995548668;
K2=   3.80023995563911;
L=   0.85120255261552;
S1=                  0;
S2=                  0;
THETA=   1.00000000000000;
V1=   0.59202988402399;
V2=   0.59202988402399;
V=   -17.6239;
end;

vcov = [ 0.001 ];

order = 6;
