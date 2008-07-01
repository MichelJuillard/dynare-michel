function oo_ = posterior_analysis(type,arg1,arg2,arg3,options_,M_,oo_)  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
    
    info = check_posterior_analysis_data(type,M_);
    SampleSize = options_.PosteriorSampleSize;
    switch info
      case 0
        disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
        error('Check the options of the estimation command...')
      case {1,2}
        MaxMegaBytes = options_.MaximumNumberOfMegaBytes;
        drsize = size_of_the_reduced_form_model(oo_.dr);
        if drsize*SampleSize>MaxMegaBytes
            drsize=0;
        end
        SampleAddress = selec_posterior_draws(SampleSize,drsize);
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
      case {4,5}
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_);
      case 6
        [ivar,vartan] = set_stationary_variables_list;
        nvar = length(ivar);
        oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan);
      otherwise
        error(['posterior_analysis:: Check_posterior_analysis_data gave a meaningless output!'])
    end
    
    
    
function oo_ = job(type,SampleSize,arg1,arg2,arg3,options_,M_,oo_,nvar,vartan)
    narg1 = 8;
    narg2 = 10;
    if ~(nargin==narg1 | nargin==narg2)
        error('posterior_analysis:: Call to function job is buggy!')
    end
    switch type
      case 'variance'
        if nargin==narg1
            [nvar,vartan,NumberOfFiles] = ...
                dsge_posterior_theoretical_covariance(SampleSize,M_,options_,oo_);
        end
        oo_ = covariance_posterior_analysis(SampleSize,M_.dname,M_.fname,...
                                            vartan,nvar,arg1,arg2,options_.mh_conf_sig,oo_);          
      case 'decomposition'
        if nargin==narg1
            [nvar,vartan,NumberOfFiles] = ...
                dsge_posterior_theoretical_variance_decomposition(SampleSize,M_,options_,oo_);
        end
        oo_ = variance_decomposition_posterior_analysis(SampleSize,M_.dname,M_.fname,...
                                                        M_.exo_names,arg2,vartan,arg1,options_.mh_conf_sig,oo_);
      case 'correlation'
        if nargin==narg1
            [nvar,vartan,NumberOfFiles] = ...
                dsge_posterior_theoretical_correlation(SampleSize,arg3,M_,options_,oo_);
        end
        oo_ = correlation_posterior_analysis(SampleSize,M_.dname,M_.fname,...
                                             vartan,nvar,arg1,arg2,arg3,options_.mh_conf_sig,oo_,M_,options_);          
      otherwise
        disp('Not yet implemented')
    end