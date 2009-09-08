var 
C
CF
CF_STAR
CH
CH_STAR
CN
CN_STAR
CT
CT_STAR
C_STAR
E
KE
KE_STAR
L
L_STAR
P
PF
PF_STAR
PH
PH_STAR
PN
PN_STAR
PT
PT_STAR
P_STAR
W
W_STAR
Y
Y_STAR
;

varexo k k_star m m_star;

parameters epsi chi thet nu phi gam;

epsi = 0.5;
nu = 3;
chi = 1.2;
phi = 4;
thet = 3;
gam = 0.5;

model;
C = (1/chi)*(exp(m)/P)^epsi;
C_STAR = (1/chi)*(exp(m_star)/P_STAR)^epsi;
CN = (1-gam)*(P/PN)*C;
CN_STAR = (1-gam)*(P_STAR/PN_STAR)*C_STAR;
CT = gam*(P/PT)*C;
CT_STAR = gam*(P_STAR/PT_STAR)*C_STAR;
CH = 0.5*(PT/PH)*CT;
CH_STAR = 0.5*(PT_STAR/PH_STAR)*CT_STAR;
CF = 0.5*(PT/PF)*CT;
CF_STAR = 0.5*(PT_STAR/PF_STAR)*CT_STAR;
P = PT^gam*PN^(1-gam);
P_STAR = PT_STAR^gam*PN_STAR^(1-gam);
PT = sqrt(PH*PF);
PT_STAR = sqrt(PH_STAR*PF_STAR);
PH = (thet/(thet-1))*W(-1);
PF_STAR = (thet/(thet-1))*W_STAR(-1);
PN = PH;
PN_STAR = PF_STAR;
L = Y;
L_STAR = Y_STAR;
(L(+1)/(P(+1)*C(+1)))*W = (phi/(phi-1))*KE(+1)*L(+1)^nu;
(L_STAR(+1)/(P_STAR(+1)*C_STAR(+1)))*W_STAR = (phi/(phi-1))*KE_STAR(+1)*L_STAR(+1)^nu;
P*C = Y*PH;
P_STAR*C_STAR = Y_STAR*PF_STAR;
Y = CH + CH_STAR + CN;
Y_STAR = CF + CF_STAR + CN_STAR;
PT = E*PT_STAR;
KE = exp(k);
KE_STAR = exp(k_star);
end;

initval;
C = 1;
PH = 1;
P = 1;
PN = 1;
PT = 1;
L = 1;
Y = 1;
W = 1;
CF = 0.25;
CH = 0.25;
CT = 0.5;
CN = 0.5;
PF = 1;
C_STAR = 1;
PH_STAR = 1;
P_STAR = 1;
PN_STAR = 1;
PT_STAR = 1;
L_STAR = 1;
Y_STAR = 1;
W_STAR = 1;
CF_STAR = 0.25;
CH_STAR = 0.25;
CT_STAR = 0.5;
CN_STAR = 0.5;
PF_STAR = 1;
KE = 1;
KE_STAR = 1;
E = 1;
k = 0;
k_star = 0;
m = 0;
m_star = 0;
end;

vcov = [
0.01 0 -0.01 0;
0 0.01 0 -0.01;
-0.01 0 0.01 0;
0 -0.01 0 0.01
];

order=4;
