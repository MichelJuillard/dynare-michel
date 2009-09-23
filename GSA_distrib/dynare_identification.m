function dynare_identification()

% main 

global M_ options_ oo_ bayestopt_ estim_params_


options_.mode_compute = 0;
[data,rawdata]=dynare_estimation_init([],1);
% computes a first linear solution to set up various variables
dynare_resolve;

options_.prior_mc=2000;

SampleSize = options_.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

prior_draw(1,bayestopt_);
IdentifDirectoryName = CheckPath('identification');

iteration = 0;
indx = estim_params_.param_vals(:,1);
indexo = estim_params_.var_exo(:,1);
useautocorr = 0;
nlags = 3;

while iteration < SampleSize,
  loop_indx = loop_indx+1;
  params = prior_draw();
  set_all_parameters(params);

  [JJ, H] = getJJ(M_,oo_,options_,0,indx,indexo,mf,nlags,useautocorr);  
  
  if ~isempty(jj),
    iteration = iteration + 1;
    pdraws(iteration,:) = params';
    [indok, indno, indweak] = identification_checks(H,JJ);    
  end
end


