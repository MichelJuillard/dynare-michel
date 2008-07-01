function oo_ = correlation_posterior_analysis(SampleSize,dname,fname,vartan,nvar,var1,var2,nar,mh_conf_sig,oo_,M_,options_)
    indx1 = check_name(vartan,var1);
    if isempty(indx1)
        disp(['posterior_analysis:: ' var1 ' is not a stationary endogenous variable!'])
        return
    end
    if ~isempty(var2)
        indx2 = check_name(vartan,var2);
        if isempty(indx2)
            disp(['posterior_analysis:: ' var2 ' is not a stationary endogenous variable!'])
            return
        end
    else
        indx2 = indx1;
        var2 = var1;
    end
    if isfield(oo_,'PosteriorTheoreticalMoments')
        if isfield(oo_.PosteriorTheoreticalMoments,'dsge')
            if isfield(oo_.PosteriorTheoreticalMoments.dsge,'correlation')
                if isfield(oo_.PosteriorTheoreticalMoments.dsge.correlation.mean,var1)
                    eval(['s1 = oo_.PosteriorTheoreticalMoments.dsge.correlation.mean' '.' var1 ';'])  
                    if isfield(s1,var2)
                        eval(['s2 = s1' '.' var2 ';'])
                        l1 = length(s2);
                        if l1<nar
                            % INITIALIZATION:
                            oo_ = initialize_output_structure(var1,var2,nar,oo_);
                            system(['rm ' M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations*']);
                            [nvar,vartan,NumberOfFiles] = ...
                                dsge_posterior_theoretical_correlation(SampleSize,nar,M_,options_,oo_);
                        else
                            if ~isnan(s2(nar))
                                %Nothing to do.
                                return
                            end
                        end
                    else
                        oo_ = initialize_output_structure(var1,var2,nar,oo_);
                    end
                else
                    oo_ = initialize_output_structure(var1,var2,nar,oo_);
                end
            else
                oo_ = initialize_output_structure(var1,var2,nar,oo_);
            end
        else
            oo_ = initialize_output_structure(var1,var2,nar,oo_);
        end
    else
        oo_ = initialize_output_structure(var1,var2,nar,oo_);
    end
    tmp = dir([ dname '/metropolis/'  fname '_PosteriorCorrelations*.mat']);
    NumberOfFiles = length(tmp);
    i1 = 1; tmp = zeros(SampleSize,1);
    for file = 1:NumberOfFiles
        load([ dname '/metropolis/'  fname '_PosteriorCorrelations' int2str(file) '.mat']);
        i2 = i1 + rows(Correlation_array) - 1;
        tmp(i1:i2) = Correlation_array(:,indx1,indx2,nar);
        i1 = i2+1;
    end
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
    name = [ var1 '.' var2 ];
    if isfield(oo_,'PosteriorTheoreticalMoments')
        if isfield(oo_.PosteriorTheoreticalMoments,'dsge')
            if isfield(oo_.PosteriorTheoreticalMoments.dsge,'correlation')
                oo_ = fill_output_structure(var1,var2,oo_,'mean',nar,post_mean);
                oo_ = fill_output_structure(var1,var2,oo_,'median',nar,post_median);
                oo_ = fill_output_structure(var1,var2,oo_,'variance',nar,post_var);
                oo_ = fill_output_structure(var1,var2,oo_,'hpdinf',nar,hpd_interval(1));
                oo_ = fill_output_structure(var1,var2,oo_,'hpdsup',nar,hpd_interval(2));
                oo_ = fill_output_structure(var1,var2,oo_,'deciles',nar,post_deciles);
                oo_ = fill_output_structure(var1,var2,oo_,'density',nar,density);
            end
        end
    end

    
function oo_ = initialize_output_structure(var1,var2,nar,oo_)
    name = [ var1 '.' var2 ];
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.mean.' name ' = NaN(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.median.' name ' = NaN(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.variance.' name ' = NaN(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.hpdinf.' name ' = NaN(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.hpdsup.' name ' = NaN(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.deciles.' name ' = cell(' int2str(nar) ',1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.density.' name ' = cell(' int2str(nar) ',1);']);
    for i=1:nar
        eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.density.' name '(' int2str(i) ',1) = {NaN};']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.deciles.' name '(' int2str(i) ',1) = {NaN};']);
    end
    
function oo_ = fill_output_structure(var1,var2,oo_,type,lag,result)
    name = [ var1 '.' var2 ];
    switch type
      case {'mean','median','variance','hpdinf','hpdsup'} 
        eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.' type '.' name '(' int2str(lag) ',1) = result;']);
      case {'deciles','density'}
        eval(['oo_.PosteriorTheoreticalMoments.dsge.correlation.' type '.' name '(' int2str(lag) ',1) = {result};']);
      otherwise
        disp('fill_output_structure:: Unknown field!')
    end