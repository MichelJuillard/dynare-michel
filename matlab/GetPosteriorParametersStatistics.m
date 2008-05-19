function get_posterior_parameters_statistics()
% This function prints and saves posterior estimates after the mcmc
% (+updates of oo_ & TeX output). 
% 
% INPUTS 
%   None.
%  
% OUTPUTS 
%   None.  
%
% SPECIAL REQUIREMENTS
%   None.
%  
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.

global estim_params_ M_ options_ bayestopt_ oo_

if ~options_.mh_replic & options_.load_mh_file
   load([M_.fname '_results.mat'],'oo_'); 
end

TeX   	= options_.TeX;
nblck 	= options_.mh_nblck;
nvx   	= estim_params_.nvx;
nvn   	= estim_params_.nvn;
ncx   	= estim_params_.ncx;
ncn   	= estim_params_.ncn;
np    	= estim_params_.np ;
nx    	= nvx+nvn+ncx+ncn+np;

DirectoryName = CheckPath('metropolis');
OutputDirectoryName = CheckPath('Output');

load([ DirectoryName '/'  M_.fname '_mh_history'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; ifil = FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2))
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
FirstMhFile = record.KeepedDraws.FirstMhFile;
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;

pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
tit2 = sprintf('%10s %7s %10s %14s %4s %6s\n',' ','prior mean','post. mean','conf. interval','prior','pstdev');
pformat = '%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f';

disp(' ');disp(' ');disp('ESTIMATION RESULTS');disp(' ');
disp(sprintf('Log data density is %f.',oo_.MarginalDensity.ModifiedHarmonicMean))
if np
    type = 'parameters';
    if TeX
        fid = TeXBegin(OutputDirectoryName,M_.fname,1,type);
    end
    disp(' ')
    disp(type)
    disp(tit2)
    ip = nvx+nvn+ncx+ncn+1;
    for i=1:np
        if options_.mh_replic
            Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, ...
             density] = posterior_moments(Draws,1,options_.mh_conf_sig);
            name = bayestopt_.name{ip};
            oo_ = Filloo(oo_,name,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
        else
            name = bayestopt_.name{ip};
            [post_mean,hpd_interval,post_var] = Extractoo(oo_,name,type);
        end
        disp(sprintf(pformat,name,bayestopt_.pmean(ip),...
                     post_mean, ...
                     hpd_interval, ...
                     pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.pstdev(ip)));    
        if TeX
            TeXCore(fid,name,deblank(pnames(bayestopt_.pshape(ip)+1,:)),bayestopt_.pmean(ip),...
                    bayestopt_.pstdev(ip),post_mean,sqrt(post_var),hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid,1,type);
    end
end
if nvx
    type = 'shocks_std';
    if TeX
        fid = TeXBegin(OutputDirectoryName,M_.fname,2,'standard deviation of structural shocks');
    end
    ip = 1;
    disp(' ')
    disp('standard deviation of shocks')
    disp(tit2)
    for i=1:nvx
        if options_.mh_replic
            Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(Draws,1,options_.mh_conf_sig);
            k = estim_params_.var_exo(i,1);
            name = deblank(M_.exo_names(k,:));
            oo_ = Filloo(oo_,name,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
            M_.Sigma_e(k,k) = post_mean*post_mean;
        else
            k = estim_params_.var_exo(i,1);
            name = deblank(M_.exo_names(k,:));
            [post_mean,hpd_interval,post_var] = Extractoo(oo_,name,type);
        end
        disp(sprintf(pformat,name,bayestopt_.pmean(ip),post_mean,hpd_interval,...
                     pnames(bayestopt_.pshape(ip)+1,:),bayestopt_.pstdev(ip)));
        if TeX
            TeXCore(fid,name,deblank(pnames(bayestopt_.pshape(ip)+1,:)),bayestopt_.pmean(ip),...
                    bayestopt_.pstdev(ip),post_mean,sqrt(post_var),hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid,2,'standard deviation of structural shocks');        
    end
end
if nvn
    type = 'measurement_errors_std';
    if TeX
        fid = TeXBegin(OutputDirectoryName,M_.fname,3,'standard deviation of measurement errors')
    end
    disp(' ')
    disp('standard deviation of measurement errors')
    disp(tit2)
    ip = nvx+1;
    for i=1:nvn
        if options_.mh_replic
            Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(Draws,1,options_.mh_conf_sig);
            name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
            oo_ = Filloo(oo_,name,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
        else
            name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
            [post_mean,hpd_interval,post_var] = Extractoo(oo_,name,type);
        end
        disp(sprintf(pformat,name,bayestopt_.pmean(ip),post_mean,hpd_interval, ...
                     pnames(bayestopt_.pshape(ip)+1,:),bayestopt_.pstdev(ip)));
        if TeX
            TeXCore(fid,name,deblank(pnames(bayestopt_.pshape(ip)+1,:)),bayestopt_.pmean(ip),...
                    bayestopt_.pstdev(ip),post_mean,sqrt(post_var),hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid,3,'standard deviation of measurement errors');        
    end
end
if ncx
    type = 'shocks_corr';
    if TeX
        fid = TeXBegin(OutputDirectoryName,M_.fname,4,'correlation of structural shocks');
    end
    disp(' ')
    disp('correlation of shocks')
    disp(tit2)
    ip = nvx+nvn+1;
    for i=1:ncx
        if options_.mh_replic
            Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(Draws,1,options_.mh_conf_sig);
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
            NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
            oo_ = Filloo(oo_,NAME,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
            M_.Sigma_e(k1,k2) = post_mean*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
            M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
        else
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
            NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
            [post_mean,hpd_interval,post_var] = Extractoo(oo_,NAME,type);
        end
        disp(sprintf(pformat, name,bayestopt_.pmean(ip),post_mean,hpd_interval, ...
                     pnames(bayestopt_.pshape(ip)+1,:),bayestopt_.pstdev(ip)));
        if TeX
            TeXCore(fid,name,deblank(pnames(bayestopt_.pshape(ip)+1,:)),bayestopt_.pmean(ip),...
                    bayestopt_.pstdev(ip),post_mean,sqrt(post_var),hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid,4,'correlation of structural shocks');
    end
end
if ncn
    type = 'measurement_errors_corr';
    if TeX
        fid = TeXBegin(OutputDirectoryName,M_.fname,5,'correlation of measurement errors');
    end
    disp(' ')
    disp('correlation of measurement errors')
    disp(tit2)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
        if options_.mh_replic
            Draws = GetAllPosteriorDraws(ip,FirstMhFile,FirstLine,TotalNumberOfMhFiles,NumberOfDraws);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(Draws,1,options_.mh_conf_sig);
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
            NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
            oo_ = Filloo(oo_,NAME,type,post_mean,hpd_interval,...
                         post_median,post_var,post_deciles,density);
        else
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
            NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
            [post_mean,hpd_interval,post_var] = Extractoo(oo_,NAME,type);
        end
        disp(sprintf(pformat, name,bayestopt_.pmean(ip),post_mean,hpd_interval, ...
                     pnames(bayestopt_.pshape(ip)+1,:),bayestopt_.pstdev(ip)));
        if TeX
            TeXCore(fid,name,deblank(pnames(bayestopt_.pshape(ip)+1,:)),bayestopt_.pmean(ip),...
                    bayestopt_.pstdev(ip),post_mean,sqrt(post_var),hpd_interval);            
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid,5,'correlation of measurement errors');        
    end
end


%
%% subfunctions:
%
function fid = TeXBegin(directory,fname,fnum,title)
    TeXfile = [directory '/' fname '_Posterior_Mean_' int2str(fnum) '.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by Dynare.\n');
    fprintf(fidTeX,['%% RESULTS FROM METROPOLIS HASTINGS (' title ')\n']);
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,['  & Prior distribution & Prior mean  & Prior ' ...
                    's.d. & Posterior mean & Posterior s.d.  & HPD inf & HPD sup\\\\ \n']);
    fprintf(fidTeX,'\\hline \\\\ \n');
    fid = fidTeX;

    
function TeXCore(fid,name,shape,priormean,priorstd,postmean,poststd,hpd)
    fprintf(fid,['$%s$ & %s & %7.3f & %6.4f & %7.3f& %6.4f & %7.4f & %7.4f \\\\ \n'],... 
            name,...
            shape,...
            priormean,...
            priorstd,...
            postmean,...
            poststd,...
            hpd(1),...
            hpd(2));
    

function TeXEnd(fid,fnum,title)
    fprintf(fid,'\\hline\\hline \n');
    fprintf(fid,'\\end{tabular}\n ');    
    fprintf(fid,['\\caption{Results from Metropolis-Hastings (' title ')}\n ']);
    fprintf(fid,['\\label{Table:MHPosterior:' int2str(fnum)  '}\n']);
    fprintf(fid,'\\end{table}\n');
    fprintf(fid,'%% End of TeX file.\n');
    fclose(fid);
    
    
function oo = Filloo(oo,name,type,postmean,hpdinterval,postmedian,postvar,postdecile,density)
    eval(['oo.posterior_mean.' type '.' name ' = postmean;']);
    eval(['oo.posterior_hpdinf.' type '.' name ' = hpdinterval(1);']); 
    eval(['oo.posterior_hpdsup.' type '.' name ' = hpdinterval(2);']);      
    eval(['oo.posterior_median.' type '.' name ' = postmedian;']);
    eval(['oo.posterior_variance.' type '.' name ' = postvar;']);
    eval(['oo.posterior_deciles.' type '.' name ' = postdecile;']);
    eval(['oo.posterior_density.' type '.' name ' = density;']);
    
function [post_mean,hpd_interval,post_var] = Extractoo(oo,name,type)
    hpd_interval = zeros(2,1);
    eval(['post_mean = oo.posterior_mean.' type '.' name ';']);
    eval(['hpd_interval(1) = oo.posterior_hpdinf.' type '.' name ';']); 
    eval(['hpd_interval(2) = oo.posterior_hpdsup.' type '.' name ';']);
    eval(['post_var = oo.posterior_variance.' type '.' name ';']);