function GetPosteriorParametersStatistics()
% stephane.adjemian@ens.fr [09-09-2005]
global estim_params_ M_ options_ bayestopt_ oo_

TeX   	= options_.TeX;
nblck 	= options_.mh_nblck;
nvx   	= estim_params_.nvx;
nvn   	= estim_params_.nvn;
ncx   	= estim_params_.ncx;
ncn   	= estim_params_.ncn;
np    	= estim_params_.np ;
nx    	= nvx+nvn+ncx+ncn+np;

DirectoryName = CheckPath('metropolis');
load([ DirectoryName '\'  M_.fname '_mh_history'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; ifil = FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2))
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
FirstMhFile = record.KeepedDraws.FirstMhFile;
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;

disp(' ')
disp(' ')
disp('ESTIMATION RESULTS')
disp(' ')
disp(sprintf('Log data density is %f.',oo_.MarginalDensity.ModifiedHarmonicMean))
pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
tit2 = sprintf('%10s %7s %10s %14s %4s %6s\n',' ','prior mean','post. mean','conf. interval','prior','pstdev');
if np
  disp(' ')
  disp('parameters')
  disp(tit2)
  ip = nvx+nvn+ncx+ncn+1;
  for i=1:np
    Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws,1);
    name = bayestopt_.name{ip};
    disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		 name, ...
		 bayestopt_.pmean(ip),post_mean,hpd_interval, ...
		 pnames(bayestopt_.pshape(ip)+1,:), ...
		 bayestopt_.pstdev(ip)));
    eval(['oo_.posterior_mean.parameters.' name ' = post_mean;']);
    eval(['oo_.posterior_hpdinf.parameters.' name ' = hpd_interval(1);']); 
    eval(['oo_.posterior_hpdsup.parameters.' name ' = hpd_interval(2);']);
    eval(['oo_.posterior_median.' name ' = post_median;']);
    eval(['oo_.posterior_variance.' name ' = post_var;']);
    eval(['oo_.posterior_deciles.' name ' = post_deciles;']);
    eval(['oo_.posterior_density.' name ' = density;']);    
    ip = ip+1;
  end
end
if nvx
  ip = 1;
  disp(' ')
  disp('standard deviation of shocks')
  disp(tit2)
  for i=1:nvx
    Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws,1);
    k = estim_params_.var_exo(i,1);
    name = deblank(M_.exo_names(k,:));
    disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		 name,bayestopt_.pmean(ip),post_mean, ...
		 hpd_interval,pnames(bayestopt_.pshape(ip)+1,:), ...
		 bayestopt_.pstdev(ip))); 
    M_.Sigma_e(k,k) = post_mean*post_mean;
    eval(['oo_.posterior_mean.shocks_std.' name ' = post_mean;']);
    eval(['oo_.posterior_hpdinf.shocks_std.' name ' = hpd_interval(1);']); 
    eval(['oo_.posterior_hpdsup.shocks_std.' name ' = hpd_interval(2);']);
    eval(['oo_.posterior_median.shocks_std.' name ' = post_median;']);
    eval(['oo_.posterior_variance.shocks_std.' name ' = post_var;']);
    eval(['oo_.posterior_deciles.shocks_std.' name ' = post_deciles;']);
    eval(['oo_.posterior_density.shocks_std.' name ' = density;']);
    ip = ip+1;
  end
end
if nvn
  disp(' ')
  disp('standard deviation of measurement errors')
  disp(tit2)
  ip = nvx+1;
  for i=1:nvn
    Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws,1);    
    name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
    disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		 name,...
		 bayestopt_.pmean(ip), ...
		 post_mean,hpd_interval, ...
		 pnames(bayestopt_.pshape(ip)+1,:), ...
		 bayestopt_.pstdev(ip)));
    eval(['oo_.posterior_mean.measurement_errors_std.' name  ' = post_mean;']);
    eval(['oo_.posterior_hpdinf.measurement_errors_std.' name ' = hpd_interval(1);']); 
    eval(['oo_.posterior_hpdsup.measurement_errors_std.' name ' = hpd_interval(2);']);		      
    eval(['oo_.posterior_median.measurement_errors_std.' name ' = post_median;']);
    eval(['oo_.posterior_variance.measurement_errors_std.' name ' = post_var;']);
    eval(['oo_.posterior_deciles.measurement_errors_std.' name ' = post_deciles;']);
    eval(['oo_.posterior_density.measurement_errors_std.' name ' = density;']);    
    ip = ip+1;
  end
end
if ncx
  disp(' ')
  disp('correlation of shocks')
  disp(tit2)
  ip = nvx+nvn+1;
  for i=1:ncx
    Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws,1);
    k1 = estim_params_.corrx(i,1);
    k2 = estim_params_.corrx(i,2);
    name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
    NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
    disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		 bayestopt_.pmean(ip),post_mean,hpd_interval, ...
		 pnames(bayestopt_.pshape(ip)+1,:), ...
		 bayestopt_.pstdev(ip)));
    eval(['oo_.posterior_mean.shocks_corr.' NAME ' = post_mean;']);
    eval(['oo_.posterior_hpdinf.shocks_corr.' NAME ' = hpd_interval(1);']); 
    eval(['oo_.posterior_hpdsup.shocks_corr.' NAME ' = hpd_interval(2);']);
    eval(['oo_.posterior_median.shocks_corr.' NAME ' = post_median;']);
    eval(['oo_.posterior_variance.shocks_corr.' NAME ' = post_var;']);
    eval(['oo_.posterior_deciles.shocks_corr.' NAME ' = post_deciles;']);
    eval(['oo_.posterior_density.shocks_corr.' NAME ' = density;']);	  
    M_.Sigma_e(k1,k2) = post_mean*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
    M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
    ip = ip+1;
  end
end
if ncn
  disp(' ')
  disp('correlation of measurement errors')
  disp(tit2)
  ip = nvx+nvn+ncx+1;
  for i=1:ncn
    Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(Draws,1);
    k1 = estim_params_.corrn(i,1);
    k2 = estim_params_.corrn(i,2);
    name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
    NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
    disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		 bayestopt_.pmean(ip),post_mean,hpd_interval, ...
		 pnames(bayestopt_.pshape(ip)+1,:), ...
		 bayestopt_.pstdev(ip))); 
    eval(['oo_.posterior_mean.measurement_errors_corr.' NAME ' = post_mean;']);
    eval(['oo_.posterior_hpdinf.measurement_errors_corr.' NAME ' = hpd_interval(1);']); 
    eval(['oo_.posterior_hpdsup.measurement_errors_corr.' NAME ' = hpd_interval(2);']);      
    eval(['oo_.posterior_median.measurement_errors_corr.' NAME ' = post_median;']);
    eval(['oo_.posterior_variance.measurement_errors_corr.' NAME ' = post_var;']);
    eval(['oo_.posterior_deciles.measurement_errors_corr.' NAME ' = post_decile;']);
    eval(['oo_.posterior_density.measurement_errors_corr.' NAME ' = density;']);
    ip = ip+1;
  end
end