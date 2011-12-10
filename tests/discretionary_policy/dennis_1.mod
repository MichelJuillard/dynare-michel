var y i pi pi_c q;
varexo g u e;

parameters omega sigma beta kappa alpha;

eta = 1;
alpha = 0.4;
beta = 0.99;
theta = 0.75;
sigma = 1;
phi = 3;

omega = 1+alpha*(2-alpha)*(sigma*eta-1);
kappa = (1-theta)*(1-beta*theta)*(phi+sigma/omega)/theta;


model(linear);
y = y(+1) -(omega/sigma)*(i-pi(+1))+g;
pi =  beta*pi(+1)+kappa*y+u;
pi_c = pi+(alpha/(1-alpha))*(q-q(-1));
q = q(+1)-(1-alpha)*(i-pi(+1))+(1-alpha)*e;
i = 1.5*pi;
end;

shocks;
var g; stderr 1;
var u; stderr 1;
var e; stderr 1;
end;

planner_objective pi_c^2 + y^2;

discretionary_policy(instruments=(i),irf=0,qz_criterium=0.999999);

