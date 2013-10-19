var y y_s R pie dq pie_s de A y_obs pie_obs R_obs vv ww;
varexo e_R e_q e_ys e_pies e_A;

parameters psi1 psi2 psi3 rho_R tau alpha rr k rho_q rho_A rho_ys rho_pies;

psi1 = 1.54;
psi2 = 0.25;
psi3 = 0.25;
rho_R = 0.5;
alpha = 0.3;
rr = 2.51;
k = 0.5;
tau = 0.5;
rho_q = 0.4;
rho_A = 0.2;
rho_ys = 0.9;
rho_pies = 0.7;

@#if !block && !bytecode
model;
@#else
@# if block && !bytecode
model(block, cutoff=0);
@# else
@#  if !block && bytecode
model(bytecode);
@#  else
model(block, bytecode, cutoff=0);
@#  endif
@# endif
@#endif
y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
pie = de+(1-alpha)*dq+pie_s;
R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;
dq = rho_q*dq(-1)+e_q;
y_s = rho_ys*y_s(-1)+e_ys;
pie_s = rho_pies*pie_s(-1)+e_pies;
A = rho_A*A(-1)+e_A;
y_obs = y-y(-1)+A;
pie_obs = 4*pie;
R_obs = 4*R;
vv = 0.2*ww+0.5*vv(-1)+1;
ww = 0.1*vv+0.5*ww(-1)+2;
/* A lt=
 0.5*vv-0.2*ww = 1
-0.1*vv+0.5*ww = 2
[ 0.5 -0.2][vv]   [1]
                =
[-0.1  0.5][ww]   [2]
det = 0.25-0.02 = 0.23
[vv]           [0.5  0.2] [1]           [0.9]   [3.91304]
     = 1/0.23*                = 1/0.23*       = 
[ww]           [0.1  0.5] [2]           [1.1]   [4.7826]
*/
end;

shocks;
var e_R = 1.25^2;
var e_q = 2.5^2;
var e_A = 1.89;
var e_ys = 1.89;
var e_pies = 1.89;
end;

options_.simul.maxit=100;
steady(solve_algo = @{solve_algo});

@#if block
model_info;
@#endif

check;

shocks;
var e_q;
periods 1;
values 0.5;
end;

simul(periods=20, markowitz=0, stack_solve_algo = @{stack_solve_algo});
/*
rplot vv;
rplot ww;
rplot A;
rplot pie;
*/

