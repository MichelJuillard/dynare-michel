function [info,description] = check_posterior_analysis_data(type,M_)
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
    
    info = 0;
    if nargout>1
        description = '';
    end
    
    %% Get informations about mcmc files. 
    if ~exist([ M_.dname '/metropolis'],'dir')
        disp('check_posterior_analysis_data:: Can''t find any mcmc file!')
        return
    end
    mhname = get_name_of_the_last_mh_file(M_);
    mhdate = get_date_of_a_file(mhname);
    
    %% Get informations about _posterior_draws files.
    if ~exist([ M_.dname '/metropolis/' M_.fname '_posterior_draws.mat'])
        info = 1; % select_posterior_draws has to be called first.
        if nargout>1
            description = 'select_posterior_draws has to be called.';
        end
        return
    else
        pddate = get_date_of_a_file([ M_.dname '/metropolis/' M_.fname '_posterior_draws.mat']);
        if pddate<mhdate
            info = 2; % _posterior_draws files have to be updated.
            if nargout>1
                description = 'posterior draws files have to be updated.';
            end
            return
        else
            info = 3; % Ok!
            if nargout>1
                description = 'posterior draws files are up to date.';
            end
        end
    end
    
    %% Get informations about posterior data files.
    switch type
      case 'variance'
        generic_post_data_file_name = 'Posterior2ndOrderMoments';
      case 'decomposition'
        generic_post_data_file_name = 'PosteriorVarianceDecomposition';
      case 'dynamic_decomposition'
        generic_post_data_file_name = 'PosteriorDynamicVarianceDecomposition';
      otherwise
        disp('This feature is not yest implemented!')
    end
    pdfinfo = dir([ M_.dname '/metropolis/' M_.fname '_' generic_post_data_file_name '*']);
    if isempty(pdfinfo)
        info = 4; % posterior draws have to be processed.
        if nargout>1
            description = 'posterior draws files have to be processed.';
        end
        return
    else
        number_of_the_last_post_data_file = length(pdfinfo);
        name_of_the_last_post_data_file = ...
            [ './' M_.dname ...
              '/metropolis/' ...
              M_.fname '_' ... 
              generic_post_data_file_name ...
              int2str(number_of_the_last_post_data_file) ...
              '.mat' ];
        pdfdate = get_date_of_a_file(name_of_the_last_post_data_file);
        if pdfdate<pddate
            info = 5; % posterior data files have to be updated.
            if nargout>1
                description = 'posterior data files have to be updated.';
            end            
        else
            info = 6; % Ok (nothing to do ;-)
            if nargout>1
                description = 'There is nothing to do';
            end        
        end
    end