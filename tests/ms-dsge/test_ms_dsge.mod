var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
phi   = 0.1;

model(use_dll);
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+3))*alpha*y(+2)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

markov_switching(chain=1, number_of_regimes=2, duration=2.5, parameters=[alpha, delta, theta], number_of_lags=1);
alpha.prior(shape=gamma, mean=3.22);
rho.prior(shape=uniform, mean=322,variance=2^.33,domain=[0.36     ,    0.88]);
std(e).prior(shape=beta,mean=0.3,variance=0.1^2,domain=[-0.1 006]);
std(y).prior(shape=beta,mean=0.3,variance=0.1^2,domain=[01 4]);
std(c).prior(shape=beta,mean=0.3,variance=0.1^2,stdev=0.2);
corr(y,c).prior(shape=beta,mean=0.3,variance=0.1^2,mode=33);
corr(b,c).prior(shape=beta,mean=0.3,variance=0.1^2);
corr(e,u).prior(shape=beta,mean=0.3,variance=0.1^2);
alpha.options(init=1);
rho.options(init=1);
beta.options(init=0.2);
std(u).options(init=3);
corr(y,c).options(init=.02);