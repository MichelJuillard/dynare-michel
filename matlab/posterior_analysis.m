function posterior_analysis(type,arg1,arg2,options_,M_,oo_)  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
    
    info = check_posterior_analysis_data(type,M_);

    switch info
      case 0
        disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
        error('Check the options of the estimation command...')
      case {1,2}
        SampleSize = options_.PosteriorSampleSize;
        MaxMegaBytes = options_.MaximumNumberOfMegaBytes;
        drsize = size_of_the_reduced_form_model(oo_.dr);
        if drsize*SampleSize>MaxMegaBytes
            drsize=0;
        end
        SampleAddress = selec_posterior_draws(SampleSize,drsize);
      case {4,5}
        switch type
          case 'variance'
            [nvar,vartan,NumberOfFiles] = ...
                dsge_posterior_theoretical_covariance(SampleSize,M_,options_,oo_);
          case 'decomposition'
            [nvar,vartan,NumberOfFiles] = ...
                dsge_posterior_theoretical_variance_decomposition(SampleSize,M_,options_,oo_);
          otherwise
            disp('Not yet implemented')
        end
      case 6
        switch type
          case 'variance'
            covariance_posterior_analysis(NumberOfFiles,SampleSize,M_.dname,M_.fname,...
                                          vartan,nvar,arg1,arg2);
          case 'decomposition'
            variance_decomposition_posterior_analysis(NumberOfFiles,SampleSize,M_.dname,M_.fname,...
                                                      M_.exo_names,arg2,vartan,nvar,arg1);
          otherwise
            disp('Not yet implemented')
        end
      otherwise
        error(['posterior_analysis:: Check_posterior_analysis_data gave a meaningless output!'])
    end
    
        
    
    
    
function covariance_posterior_analysis(NumberOfFiles,NumberOfSimulations,dname,fname,vartan,nvar,var1,var2)
    indx1 = check_name(vartan,var1)
    if isempty(indx1)
        disp(['posterior_analysis:: ' var1 ' is not a stationary endogenous variable!'])
        return
    end
    if ~isempty(var2)
        indx2 = check_name(vartan,var2)
        if isempty(indx2)
            disp(['posterior_analysis:: ' var2 ' is not a stationary endogenous variable!'])
            return
        end
    else
        indx2 = indx1;
        var2 = var1;
    end
    i1 = 1; tmp = zeros(NumberOfSimulations,1);
    for file = 1:CovarFileNumber
        load([ dname '/metropolis/'  fname '_Posterior2ndOrderMoments' int2str(file)]);
        i2 = i1 + rows(Covariance_matrix) - 1;
        tmp(i1:i2) = Covariance_matrix(:,symmetric_matrix_index(indx1,indx2,nvar));
        i1 = i2+1;
    end
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
        posterior_moments(tmp,1,options_.mh_conf_sig);
    if strcmpi(var1,var2)
        name = var1;
    else
        name = [var1 '.' var2];
    end
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.mean.' name ' = post_mean;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.median.' name ' = post_median;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.variance.' name ' = post_var;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdinf.' name ' = hpd_interval(1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdsup.' name ' = hpd_interval(2);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.deciles.' name ' = post_deciles;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.density.' name ' = density;']);
    
    
function variance_decomposition_posterior_analysis(NumberOfFiles,NumberOfSimulations,dname,fname, ...
                                          exonames,exo,vartan,nvar,var)
    indx1 = check_name(vartan,var)
    if isempty(indx)
        disp(['posterior_analysis:: ' var ' is not a stationary endogenous variable!'])
        return
    end
    jndx = check_name(exonames,exo);
    if isempty(jndx)
        disp(['posterior_analysis:: ' exo ' is not a declared exogenous variable!'])
        return
    end
    i1 = 1; tmp = zeros(NumberOfSimulations,1);
    for file = 1:NumberOfFiles
        load([fname '_PosteriorVarianceDecomposition' int2str(file)]);
        i2 = i1 + rows(Decomposition_array) - 1;
        tmp(i1:i2) = Decomposition_array(:,nvar+(i-1)*nexo+j);
        i1 = i2+1;
    end
    name = [ var '.' exo ];
        t1 = min(tmp); t2 = max(tmp);
        t3 = t2-t1;% How to normalize ? t1 and t2 may be zero...
        if t3<1.0e-12
            if t1<1.0e-12
                t1 = 0;
            end
            if abs(t1-1)<1.0e-12
                t1 = 1;
            end 
            post_mean = t1;
            post_median = t1;
            post_var = 0;
            hpd_interval = NaN(2,1);
            post_deciles = NaN(9,1);
            density = NaN;
        else
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(tmp,1,options_.mh_conf_sig);
        end
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean.' name ' = post_mean;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.median.' name ' = post_median;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.variance.' name ' = post_var;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdinf.' name ' = hpd_interval(1);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdsup.' name ' = hpd_interval(2);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.deciles.' name ' = post_deciles;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.density.' name ' = density;']);

    
function n = check_name(vartan,varname)
    n = strmatch(varname,vartan,'exact')