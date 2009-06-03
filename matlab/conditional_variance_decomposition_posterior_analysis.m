function oo_ = conditional_variance_decomposition_posterior_analysis(NumberOfSimulations, dname, fname, ...
                                          Steps, exonames, exo, vartan, var, mh_conf_sig, oo_)
% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

    indx = check_name(vartan,var);
    if isempty(indx)
        disp(['posterior_analysis:: ' var ' is not a stationary endogenous variable!'])
        return
    end
    endogenous_variable_index = sum(1:indx);
    exogenous_variable_index = check_name(exonames,exo);
    if isempty(exogenous_variable_index)
        disp(['posterior_analysis:: ' exo ' is not a declared exogenous variable!'])
        return
    end
    
    tmp = dir([ dname '/metropolis/'  fname '_PosteriorConditionalVarianceDecomposition*.mat']);
    NumberOfFiles = length(tmp);
    i1 = 1; tmp = zeros(NumberOfSimulations,length(Steps));
    
    for file = 1:NumberOfFiles
        load([dname '/metropolis/' fname '_PosteriorConditionalVarianceDecomposition' int2str(file) '.mat']);
        % (endovar,time,exovar,simul)
        i2 = i1 + size(Conditional_decomposition_array,4) - 1;
        tmp(i1:i2,:) = transpose(dynare_squeeze(Conditional_decomposition_array(endogenous_variable_index,:,exogenous_variable_index,:)));
        i1 = i2+1;
    end
    name = [ var '.' exo ];
    if isfield(oo_,'PosteriorTheoreticalMoments')
        if isfield(oo_.PosteriorTheoreticalMoments,'dsge')
            if isfield(oo_.PosteriorTheoreticalMoments.dsge,'ConditionalVarianceDecomposition')
                if isfield(oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean,name)
                    if sum(Steps-oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean.(name)(1,:)) == 0
                        % Nothing (new) to do here...
                        return
                    end
                end
            end
        end
    end
    posterior_mean = NaN(2,length(Steps));
    posterior_mean(1,:) = Steps;
    posterior_median = NaN(1,length(Steps));
    posterior_variance = NaN(1,length(Steps));
    posterior_deciles = NaN(9,length(Steps));
    posterior_density = NaN(2^9,2,length(Steps));
    posterior_hpdinf = NaN(1,length(Steps));
    posterior_hpdsup = NaN(1,length(Steps));
    for i=1:length(Steps)
        if ~isconst(tmp(:,i))
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
                posterior_moments(tmp(:,i),1,mh_conf_sig);
            posterior_mean(2,i) = post_mean;
            posterior_median(i) = post_median;
            posterior_variance(i) = post_var;
            posterior_deciles(:,i) = post_deciles;
            posterior_hpdinf(i) = hpd_interval(1);
            posterior_hpdinf(i) = hpd_interval(2);
            posterior_density(:,:,i) = density;
        end
    end
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.mean.' name ' = posterior_mean;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.median.' name ' = posterior_median;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.variance.' name ' = posterior_variance;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.hpdinf.' name ' = posterior_hpdinf;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.hpdsup.' name ' = posterior_hpdsup;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.deciles.' name ' = posterior_deciles;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.density.' name ' = posterior_density;']);