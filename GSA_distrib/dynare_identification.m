function [pdraws, idemodel, idemoments] = dynare_identification()

% main 

global M_ options_ oo_ bayestopt_ estim_params_


options_ = set_default_option(options_,'datafile',[]);
options_.mode_compute = 0;
[data,rawdata]=dynare_estimation_init([],1);
% computes a first linear solution to set up various variables
dynare_resolve;

options_.prior_mc=2000;

SampleSize = options_.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

prior_draw(1,bayestopt_);
IdentifDirectoryName = CheckPath('identification');

indx = estim_params_.param_vals(:,1);
indexo=[];
if ~isempty(estim_params_.var_exo)
  indexo = estim_params_.var_exo(:,1);
end
useautocorr = 0;
nlags = 3;

iteration = 0;
loop_indx = 0;

h = waitbar(0,'Monte Carlo identification checks ...');

while iteration < SampleSize,
  loop_indx = loop_indx+1;
  params = prior_draw();
  set_all_parameters(params);

  [JJ, H] = getJJ(M_,oo_,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);  
  
  if ~isempty(JJ),
    iteration = iteration + 1;
    pdraws(iteration,:) = params';
    [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), ...
      idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), ...
      idemodel.cond(iteration), idemoments.cond(iteration), ...
      idemodel.ee(:,iteration), idemoments.ee(:,iteration), ...
      idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
      idemodel.indno{iteration}, idemoments.indno{iteration}] = ...
      identification_checks(H,JJ, bayestopt_);   
    
      waitbar(iteration/SampleSize,h)
  end
end

close(h)

      save([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments')



