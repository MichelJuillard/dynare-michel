var y i pi rbar ;

varexo r tauw taus taua gn;

parameters khia khiw khis phipi phiy taubs taubw tauba w sigma psi kappa alpha mu beta teta;

teta = 12.7721;
sigma = 1.1599;
beta = 0.9970;
alpha = 0.7747;
mu = 0.9030;
taubs = 0.05;
taubw = 0.02;
tauba = 0;
w = 1.5692;
phipi = 1.5;
phiy = 0.5/4;

khia = (1-beta)/(1-tauba);
khiw = 1/(1-taubw);
khis = 1/(1+taubs);
psi = 1/(sigma + w);
kappa = (1-alpha)*(1-alpha*beta)*(sigma+w)/(alpha*(1+w*teta));

model(linear);
y = y(+1)-sigma*(i-pi(+1)-r)+(gn-gn(+1))+(sigma)^-1*khis*(taus(+1)-taus)+sigma*khia*taua;

pi=kappa*y+kappa*psi*(khiw*tauw+khis*taus-sigma*gn)+beta*pi(+1);

i=max(0,r+phipi*pi+phiy*y);

rbar = -((kappa*phipi+(1-beta*mu)*phiy)*sigma^-1*khia*taus)/((1-mu+sigma^-1*phiy)*(1-beta*mu)+kappa*sigma^-1*(phipi-mu))
- (((1-mu)*kappa*psi*phipi+sigma^-1*mu*kappa*psi*phiy)*khiw*tauw)/((1-mu+sigma^-1*phiy)*(1-beta*mu)+kappa*sigma*(phipi-mu))
-(kappa*sigma*(1-mu)*(sigma^-1-psi)*phipi+((1-mu)*(1-beta*mu)-kappa*psi*mu)*phiy)*(gn-sigma^-1*khis*taus)/((1-mu-sigma^-1*phiy)*(1-beta*mu)+kappa*sigma^-1*(phipi-mu));

end;

initval;
y=0;
i=-log(beta);
pi=0;
rbar = 0;
end;

steady;
check;

shocks;
var r;
periods 1:9;
values -0.0104;
end;

simul(periods=2100);
