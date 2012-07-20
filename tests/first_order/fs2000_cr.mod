// See fs2000.mod in the examples/ directory for details on the model

@#define countries = 1:100

var 
@#for c in countries
 m_@{c} P_@{c} c_@{c} e_@{c} W_@{c} R_@{c} k_@{c} d_@{c} n_@{c} l_@{c} gy_obs_@{c} gp_obs_@{c} y_@{c} dA_@{c}
@#endfor
;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
@#for c in countries
dA_@{c} = exp(gam+e_a);
log(m_@{c}) = (1-rho)*log(mst) + rho*log(m_@{c}(-1))+e_m;
-P_@{c}/(c_@{c}(+1)*P_@{c}(+1)*m_@{c})+bet*P_@{c}(+1)*(alp*exp(-alp*(gam+log(e_@{c}(+1))))*k_@{c}^(alp-1)*n_@{c}(+1)^(1-alp)+(1-del)*exp(-(gam+log(e_@{c}(+1)))))/(c_@{c}(+2)*P_@{c}(+2)*m_@{c}(+1))=0;
W_@{c} = l_@{c}/n_@{c};
-(psi/(1-psi))*(c_@{c}*P_@{c}/(1-n_@{c}))+l_@{c}/n_@{c} = 0;
R_@{c} = P_@{c}*(1-alp)*exp(-alp*(gam+e_a))*k_@{c}(-1)^alp*n_@{c}^(-alp)/W_@{c};
1/(c_@{c}*P_@{c})-bet*P_@{c}*(1-alp)*exp(-alp*(gam+e_a))*k_@{c}(-1)^alp*n_@{c}^(1-alp)/(m_@{c}*l_@{c}*c_@{c}(+1)*P_@{c}(+1)) = 0;
c_@{c}+k_@{c} = exp(-alp*(gam+e_a))*k_@{c}(-1)^alp*n_@{c}^(1-alp)+(1-del)*exp(-(gam+e_a))*k_@{c}(-1);
P_@{c}*c_@{c} = m_@{c};
m_@{c}-1+d_@{c} = l_@{c};
e_@{c} = exp(e_a);
y_@{c} = k_@{c}(-1)^alp*n_@{c}^(1-alp)*exp(-alp*(gam+e_a));
gy_obs_@{c} = dA_@{c}*y_@{c}/y_@{c}(-1);
gp_obs_@{c} = (P_@{c}/P_@{c}(-1))*m_@{c}(-1)/dA_@{c};
@#endfor
end;

initval;
@#for c in countries
k_@{c} = 6;
m_@{c} = mst;
P_@{c} = 2.25;
c_@{c} = 0.45;
e_@{c} = 1;
W_@{c} = 4;
R_@{c} = 1.02;
d_@{c} = 0.85;
n_@{c} = 0.19;
l_@{c} = 0.86;
y_@{c} = 0.6;
gy_obs_@{c} = exp(gam);
gp_obs_@{c} = exp(-gam);
dA_@{c} = exp(gam);
@#endfor
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

tic;
check;
disp(toc);

tic;
stoch_simul(order=1,dr=cycle_reduction,irf=0,nomoments,noprint);
disp(toc);
