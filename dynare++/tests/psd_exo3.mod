var c k x;
varexo rho;

parameters a alph gam bet lamb;
alph = 0.7;
bet = 0.95;
gam = 2;
a = 1.052632;
lamb = 0.9;

model;
c^(-gam) = bet*c(+1)^(-gam)*a*exp(x(+1))*k^(-alph);
k = a*exp(x)*k(-1)^(1-alph)/(1-alph)-c;
x = lamb*x(-1)+rho;
end;

initval;
k = 1;
c = 2.508;
x = 0;
rho = 0;
end;

vcov=[0.0001];

order=6;