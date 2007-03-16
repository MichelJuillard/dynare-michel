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
opt_gsa.Nsam=2048; % sample size (default = 2048)
opt_gsa.redform=1; % prepares mapping of reduced form coefficients (default = 0): this saves the full MC sample of the reduced form LRE solution
opt_gsa.ilptau=1; % lptau sample (1=default)
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
opt_gsa.namendo=['pie'; 'R  ']; % evaluate relationships for pie and R (it can be M_.endo_names as well for complete analysis)
opt_gsa.namexo='e_R'; % evaluate relationships with exogenous e_R
//opt_gsa.namexo=M_.exo_names; % evaluate relationships with all exogenous
opt_gsa.namlagendo='R'; % evaluate relationships with lagged endogenous
//opt_gsa.namlagendo=M_.endo_names; % evaluate relationships with all lagged endogenous
opt_gsa.load_stab=1; % load stability analsis sample
opt_gsa.load_redform=0; %load reduced form analysis if available (default=0: perform a new one)
opt_gsa.stab=0; % don't do again stability analysis
options_.opt_gsa=opt_gsa;
dynare_sensitivity;

disp(' ');
disp('I ESTIMATE THE MODEL');
disp(' ');
pause;
// if already estimated, use this to build filtered variables at the mode in oo_ for RMSE analysis
estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2, mode_file=ls2003_mode, load_mh_file,
  prefilter=1,mh_jscale=0.55,mh_replic=0, mode_compute=0, nograph, mh_drop=0.6);

// run this to generate posterior mode and Metropolis files if not yet done
//estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2,prefilter=1,mh_jscale=0.5,mh_replic=100000, mode_compute=4, nograph, mh_drop=0.6);
//options_.hess=1;
//options_.ftol=1.e-7;


// run this to produce posterior samples of filtered, smoothed amd irf variables
//estimation(datafile=data_ca1,first_obs=8,nobs=79,mh_nblocks=2,prefilter=1,mh_jscale=0.5,
//          mh_replic=0, mode_file=ls2003_mode, mode_compute=0, nograph, load_mh_file, bayesian_irf,
//		  filtered_vars, smoother, mh_drop=0.6);


disp(' ');
disp('MC FILTERING(opt_gsa.rmse=1), TO MAP THE FIT FROM PRIORS');
pause;

opt_gsa.redform=0;
opt_gsa.load_stab=1; % load prior sample
opt_gsa.load_rmse=0; % make a new rmse analysis
opt_gsa.istart_rmse=2; %start computing rmse from second observation (i.e. rmse does not inlude initial big error)
opt_gsa.stab=0; % don't  plot again stability analysis results
opt_gsa.rmse=1; % do rmse analysis
opt_gsa.glue=1; % prepare for glue GUI
opt_gsa.pfilt_rmse=0.1;  % filtering criterion, i.e. I filter the best 10% rmse's
opt_gsa.alpha2_rmse=0.3; % critical value for correlations in the rmse filterting analysis: if ==1, means no corrleation analysis done
opt_gsa.alpha_rmse=1; % critical value for smirnov statistics
options_.opt_gsa=opt_gsa;
dynare_sensitivity;

disp(' ');
disp('WE DO STABILITY MAPPING AGAIN, BUT FROM THE POSTERIOR RANGES (opt_gsa.pprior=0 & opt_gsa.ppost=0)');
pause,

opt_gsa.stab=1;
opt_gsa.Nsam=2048;
opt_gsa.redform=0;
opt_gsa.pprior=0;
opt_gsa.rmse=0;
opt_gsa.load_stab=0;
opt_gsa.alpha2_stab=0.4;
options_.opt_gsa=opt_gsa;
dynare_sensitivity;
//stab_map_(2048,1,0.4,1,0);

disp(' ');
disp('RMSE ANALYSIS FOR POSTERIOR RANGES');
pause,
opt_gsa.stab=0;
opt_gsa.redform=0;
opt_gsa.rmse=1;
opt_gsa.load_rmse=0;
opt_gsa.alpha2_rmse=1;
opt_gsa.alpha_rmse=1;
options_.opt_gsa=opt_gsa;
dynare_sensitivity;

disp(' ');
disp('RMSE ANALYSIS FOR POSTERIOR MCMC sample (opt_gsa.ppost=1)');
pause,
opt_gsa.stab=0;
opt_gsa.redform=0;
opt_gsa.rmse=1;
opt_gsa.load_rmse=0;
opt_gsa.ppost=1;
opt_gsa.alpha2_rmse=1;
opt_gsa.alpha_rmse=1;
opt_gsa.glue=1;
options_.opt_gsa=opt_gsa;
dynare_sensitivity;


