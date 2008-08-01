function oo_ = covariance_posterior_analysis(NumberOfSimulations,dname,fname,vartan,nvar,var1,var2,mh_conf_sig,oo_)

% Copyright (C) 2008 Dynare Team
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
            if isfield(oo_.PosteriorTheoreticalMoments.dsge,'covariance')
                if isfield(oo_.PosteriorTheoreticalMoments.dsge.covariance.mean,var1)
                    eval(['s1 = oo_.PosteriorTheoreticalMoments.dsge.covariance.mean' '.' var1 ';'])  
                    if isfield(s1,var2)
                        % Nothing to do.
                        return
                    end
                else
                    if isfield(oo_.PosteriorTheoreticalMoments.dsge.covariance.mean,var2)
                        eval(['s2 = oo_.PosteriorTheoreticalMoments.dsge.covariance.mean' '.' var2 ';'])
                        if isfield(s1,var1)
                            % Nothing to do (the covariance matrix is symmetric!).
                            return
                        end
                    end
                end         
            end
        end
    end
    tmp = dir([ dname '/metropolis/'  fname '_Posterior2ndOrderMoments*.mat']);
    NumberOfFiles= length(tmp);
    i1 = 1; tmp = zeros(NumberOfSimulations,1);
    for file = 1:NumberOfFiles
        load([ dname '/metropolis/'  fname '_Posterior2ndOrderMoments' int2str(file) '.mat']);
        i2 = i1 + rows(Covariance_matrix) - 1;
        tmp(i1:i2) = Covariance_matrix(:,symmetric_matrix_index(indx1,indx2,nvar));
        i1 = i2+1;
    end
    [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
    name = [var1 '.' var2];
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.mean.' name ' = post_mean;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.median.' name ' = post_median;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.variance.' name ' = post_var;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdinf.' name ' = hpd_interval(1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdsup.' name ' = hpd_interval(2);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.deciles.' name ' = post_deciles;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.density.' name ' = density;']);