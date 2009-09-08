var C K r w N tau I;
varexo e;

parameters alph bet delt thet tau_m rho;
alph = 0.3;
bet = 0.96;
thet = 0.3;
delt = 0.05;
tau_m = 0.35;
rho = 0.8;

model;
C = C(+1)/(bet*(r(+1)+1-delt));
I = K(-1)^alph*N^(1-alph)-C;
K = I+(1-delt)*K(-1);
N = 1-(1-thet)*C/(thet*w);
r = (1-tau)*alph*(K(-1)/N)^(alph-1);
w = (1-tau)*(1-alph)*(K(-1)/N)^alph;
tau = (1-rho)*tau_m + rho*tau(-1)+e;
end;

initval;
C=0.2;
I=0.02;
K=0.5;
N=0.18;
r=0.09;
w=0.6;
tau=0.35;
e=0;
end;

vcov = [0.007208];

order=7;





