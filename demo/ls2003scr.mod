var y y_s R pie dq pie_s de A y_obs pie_obs R_obs;
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


model(linear);
y = y(+1) - (tau +alpha*(2-alpha)*(1-tau))*(R-pie(+1))-alpha*(tau +alpha*(2-alpha)*(1-tau))*dq(+1) + alpha*(2-alpha)*((1-tau)/tau)*(y_s-y_s(+1))-A(+1);
pie = exp(-rr/400)*pie(+1)+alpha*exp(-rr/400)*dq(+1)-alpha*dq+(k/(tau+alpha*(2-alpha)*(1-tau)))*y+k*alpha*(2-alpha)*(1-tau)/(tau*(tau+alpha*(2-alpha)*(1-tau)))*y_s;
pie = de+(1-alpha)*dq+pie_s;
R = rho_R*R(-1)+(1-rho_R)*(psi1*pie+psi2*(y+alpha*(2-alpha)*((1-tau)/tau)*y_s)+psi3*de)+e_R;
dq = rho_q*dq(-1)+e_q;
y_s = rho_ys*y_s(-1)+e_ys;
pie_s = rho_pies*pie_s(-1)+e_pies;
A = rho_A*A(-1)+e_A;
y_obs = y-y(-1)+A;
pie_obs = 4*pie;
R_obs = 4*R;
end;

shocks;
	var e_R = 1.25^2;
	var e_q = 2.5^2;
	var e_A = 1.89;
	var e_ys = 1.89;
	var e_pies = 1.89;
end;

varobs y_obs R_obs pie_obs dq de;

//addpath H:\Junior\2006\gsautilities\GSA;

estimated_params;
psi1 , gamma_pdf,1.5,0.5;
psi2 , gamma_pdf,0.25,0.125;
psi3 , gamma_pdf,0.25,0.125;
rho_R ,beta_pdf,0.5,0.2;
alpha ,beta_pdf,0.3,0.1;
rr ,gamma_pdf,2.5,1;
k , gamma_pdf,0.5,0.25;
tau ,gamma_pdf,0.5,0.2;
rho_q ,beta_pdf,0.4,0.2;
rho_A ,beta_pdf,0.5,0.2;
rho_ys ,beta_pdf,0.8,0.1;
rho_pies,beta_pdf,0.7,0.15;
stderr e_R,inv_gamma_pdf,1.2533,0.6551;
stderr e_q,inv_gamma_pdf,2.5066,1.3103;
stderr e_A,inv_gamma_pdf,1.2533,0.6551;
stderr e_ys,inv_gamma_pdf,1.2533,0.6551;
stderr e_pies,inv_gamma_pdf,1.88,0.9827;
end;

  
disp(' ');
disp('NOW I DO STABILITY MAPPING, WHICH REQUIRES dynare_estimation to initialise prior settings');
disp(' ');
pause;

estimation(datafile=data_ca1,mode_compute=0);

		   
opt_gsa.stab=1;	% performs stability analysis (1 = default)
opt_gsa.morris = 1; % performs screening (default = 0)
opt_gsa.morris_nliv = 6; % number of bins in each prior range (default = 6)
opt_gsa.morris_ntra = 20; % number of replicas (default = 20)
opt_gsa.redform=1; % prepares mapping of reduced form coefficients (default = 0): this saves the full MC sample of the reduced form LRE solution
opt_gsa.load_stab=0; % generate a new sample: overwrites any generated sample (default=0)
opt_gsa.alpha2_stab=0.4; % critical value to plot correlations in stable samples (default = 0.3)
opt_gsa.ksstat=0; % critical value to plot Smirnov test in filtered samples (default = 0.1)

options_.opt_gsa=opt_gsa;
dynare_sensitivity;

disp(' ');
disp('ANALYSIS OF REDUCED FORM COEFFICIENTS');
disp(' ');
pause;

opt_gsa.logtrans_redform=1; % also estimate log-transformed reduced form coefficients (default=0)
//opt_gsa.namendo='pie'; % evaluate relationships for pie (it can be M_.endo_names as well for complete analysis)
//opt_gsa.namendo=['pie'; 'R  ']; % evaluate relationships for pie and R (it can be M_.endo_names as well for complete analysis)
opt_gsa.namendo=M_.endo_names; % evaluate relationships with all lagged endogenous
//opt_gsa.namexo='e_R'; % evaluate relationships with exogenous e_R
opt_gsa.namexo=M_.exo_names; % evaluate relationships with all exogenous
//opt_gsa.namlagendo='R'; % evaluate relationships with lagged endogenous
opt_gsa.namlagendo=M_.endo_names; % evaluate relationships with all lagged endogenous
opt_gsa.load_stab=1; % load stability analsis sample
opt_gsa.stab=0; % don't do again stability analysis
options_.opt_gsa=opt_gsa;
dynare_sensitivity;

