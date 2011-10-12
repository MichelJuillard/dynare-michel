var y c k m n R pi z u;
varexo	e sigma;	 
// sigma stands for phi in the eq 2.37 p.69

parameters alpha beta delta gamm phi1 eta a b rho  phi2 Psi thetass;  
//phi1 stands for capital phi in eq.2.68 and 2.69
//phi2 stands for lowercase phi in eq. 2.66

alpha = 0.36;
beta = 0.989; 
gamm = 0.5;
delta = 0.019;
phi1 = 2;
phi2 = 0;
eta = 1;
a = 0.95;
b = 2.56;
rho = 0.95;
Psi = 1.47630583;
thetass = 1.0125;

model;

(a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*a*exp(c)^(-b) = (a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*(1-a)*exp(m)^(-b)+beta*(a*exp(c(+1))^(1-b)+(1-a)*exp(m(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c(+1))^(-b)/(1+pi(+1));

Psi*(1-exp(n))^(-eta)/(a*exp(c)^(-b)*(a*exp(c)^(1-b) + (1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))) = (1-alpha)*exp(y)/exp(n);

(a*exp(c)^(1-b)+(1-a)*exp(m)^(1-b))^((b-phi1)/(1-b))*a*exp(c)^(-b) = beta*exp(R(+1))*(a*exp(c(+1))^(1-b)+(1-a)*exp(m(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c(+1))^(-b);

exp(R) = alpha*exp(y)/exp(k(-1)) + 1-delta;

exp(k) = (1-delta)*exp(k(-1))+exp(y)-exp(c);

exp(y) = exp(z)*exp(k(-1))^alpha*exp(n)^(1-alpha);

exp(m) = exp(m(-1))*(u+thetass)/(1+pi);

z = rho*z(-1) + e;

u = gamm*u(-1) + phi2*z(-1) + sigma;

end;

shocks;
var e; stderr 0.007;
var sigma;stderr 0.0089;
end;

steady_state_model;
// solving in levels
// calibrating n = 1/3 and recovering the value of Psi
// adapting solution Walsh (2003) p. 84
pi = thetass-1;
en = 1/3;
eR = 1/beta;
y_k = (1/alpha)*(1/beta-1+delta);
ek = en*y_k^(-1/(1-alpha));
ec = ek*(y_k-delta);
em = ec*(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^(-1/b);
ey = ek*y_k;
Xss = a*ec^(1-b)*(1+(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^((b-1)/b));
Psi = (1-alpha)*(ey/en)*Xss^((b-phi1)/(1-b))*a*ec^(-b)*(1-en)^eta;
n = log(en);
k = log(ek);
m = log(em);
c = log(ec);
y = log(ey);
R = log(eR);
end;

steady;

