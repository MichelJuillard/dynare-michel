var y c k i l y_l w r z;
varexo e;
parameters beta psi delta alpha rho epsilon;

model;
(1/c) = beta*(1/c(+1))*(1+r(+1)-delta);
psi*c/(1-l) = w;
c+i = y;
y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha);
w = y*((epsilon-1)/epsilon)*(1-alpha)/l;
r = y*((epsilon-1)/epsilon)*alpha/k(-1);
i = k-(1-delta)*k(-1);
y_l = y/l;
z = rho*z(-1)+e;
end;

varobs y;

initval;
  k = 9;
  c = 0.76;
  l = 0.3;
  w = 2.07;
  r = 0.03;
  z = 0;
  e = 0;
end;

estimated_params;
alpha, beta_pdf, 0.35, 0.02;
beta, beta_pdf, 0.99, 0.002;
delta, beta_pdf, 0.025, 0.003;
psi, gamma_pdf, 1.75, 0.02;
rho, beta_pdf, 0.95, 0.05;
epsilon, gamma_pdf, 10, 0.003;
stderr e, inv_gamma_pdf, 0.01, inf;
end;

estimation(datafile=simuldataRBC,nobs=200,first_obs=500,mh_replic=2000,mh_nblocks=2,mh_drop=0.45,mh_jscale=0.8,mode_compute=6);
