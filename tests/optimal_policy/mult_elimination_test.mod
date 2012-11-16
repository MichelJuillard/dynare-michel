// Test of mult_elimination function

// parameters value is set to posterior mean as computed by ls1.mod
// dR in objective function

var A de dq dR pie pie_obs pie_s R R_obs y y_obs y_s ;
varexo e_A e_pies e_q e_ys ;

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

model(linear);
y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
pie = de+(1-alpha)*dq+pie_s;
//R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;
dq = rho_q*dq(-1)+e_q;
y_s = rho_ys*y_s(-1)+e_ys;
pie_s = rho_pies*pie_s(-1)+e_pies;
A = rho_A*A(-1)+e_A;
y_obs = y-y(-1)+A;
pie_obs = 4*pie;
R_obs = 4*R;
dR = R-R(-1);
end;

shocks;
var e_q = 2.5^2;
var e_A = 1.89;
var e_ys = 1.89;
var e_pies = 1.89;
end;

planner_objective 0.25*pie_obs^2+y^2+0.1*dR^2;

ramsey_policy(order=1,irf=0,planner_discount=0.95);

dr2 = mult_elimination({'R'},M_,options_,oo_);

k1 = M_.nstatic+(1:M_.nspred);
k2 = strmatch('MULT_',M_.endo_names(oo_.dr.order_var(k1),:));
k3 = k1(setdiff(1:M_.nspred,k2));
k4 = oo_.dr.order_var(k3);

V0 = oo_.var(k4,k4);


Atest = [dr2.M1(k3,:) dr2.M2(k3,:) dr2.M4(k3,:); eye(6) zeros(6,10);zeros(4,16)];
Btest = [dr2.M3(k3,:); zeros(6,4); eye(4)];

V1=lyapunov_symm(Atest,Btest*M_.Sigma_e*Btest',1+1e-6,options_.lyapunov_complex_threshold);

if max(max(abs(V1(1:6,1:6)-V0)))
   disp('Test OK')
end