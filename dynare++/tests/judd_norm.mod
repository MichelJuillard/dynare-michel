var C K1 K2 L S1 S2 THETA V V1 V2;

varexo KSI;

parameters theta_ss lambda delta1 delta2 alpha1 alpha2 eta beta gamma depr1 depr2;

theta_ss=1;
lambda=0.5;
delta1=0.05;
delta2=0.2;
alpha1=0.3;
alpha2=0.3;
eta=3;
beta=0.95;
gamma=0.5;
depr1=0.1;
depr2=0.05;

model;
1 = (THETA*K1^alpha1*K2^alpha2*L^(1-alpha1-alpha2)-S1*K1-S2*K2)/C;
1 = (1-depr1+(1-0.5*delta1*S1)*S1)*K1(-1)/K1;
1 = (1-depr2+(1-0.5*delta2*S2)*S2)*K2(-1)/K2;
1 = THETA(-1)^lambda/THETA*theta_ss^(1-lambda)*exp(KSI);
/*
1 = (THETA(-1)*lambda+theta_ss*(1-lambda)+KSI)/THETA;
*/
C^(-gamma)*THETA*K1^alpha1*K2^alpha2*L^(-alpha1-alpha2)*(1-alpha1-alpha2)*L^(-eta)=1;
1 = beta*V1(+1)*(1-delta1*S1)*C^gamma; 
1 = beta*V2(+1)*(1-delta2*S2)*C^gamma;
1 = (C^(-gamma)*(alpha1*THETA*K1^(alpha1-1)*K2^alpha2*L^(1-alpha1-alpha2)-S1)+beta*V1(+1)*(1-depr1+(1-0.5*delta1*S1)*S1))/V1; 
1 = (C^(-gamma)*(alpha2*THETA*K1^alpha1*K2^(alpha2-1)*L^(1-alpha1-alpha2)-S2)+beta*V2(+1)*(1-depr2+(1-0.5*delta2*S2)*S2))/V2;
1 = (C^(1-gamma)/(1-gamma)-L^(1+eta)/(1+eta) + beta*V(+1))/V;
end;

initval;
C      =1.0997055 ;
L      =0.9425540 ;
S1     =0.1005051 ;
S2     =0.0500627 ;
K1     =2.9378521 ;
K2     =2.1952681 ;
THETA  =1.        ;
V      =38.000392 ;
V1     =1.0139701 ;
V2     =1.0062981 ;
end;

vcov = [ 0.05 ];

order = 5;
