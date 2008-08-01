function oo_ = variance_decomposition_posterior_analysis(NumberOfSimulations,dname,fname, ...
                                          exonames,exo,vartan,var,mh_conf_sig,oo_)

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

    indx = check_name(vartan,var);
    if isempty(indx)
        disp(['posterior_analysis:: ' var ' is not a stationary endogenous variable!'])
        return
    end
    jndx = check_name(exonames,exo);
    if isempty(jndx)
        disp(['posterior_analysis:: ' exo ' is not a declared exogenous variable!'])
        return
    end
    tmp = dir([ dname '/metropolis/'  fname '_PosteriorVarianceDecomposition*.mat']);
    NumberOfFiles = length(tmp);
    i1 = 1; tmp = zeros(NumberOfSimulations,1);
    indice = (indx-1)*rows(exonames)+jndx;
    for file = 1:NumberOfFiles
        load([dname '/metropolis/' fname '_PosteriorVarianceDecomposition' int2str(file) '.mat']);
        i2 = i1 + rows(Decomposition_array) - 1;
        tmp(i1:i2) = Decomposition_array(:,indice);
        i1 = i2+1;
    end
    name = [ var '.' exo ];
    if isfield(oo_,'PosteriorTheoreticalMoments')
        if isfield(oo_.PosteriorTheoreticalMoments,'dsge')
            if isfield(oo_.PosteriorTheoreticalMoments.dsge,'VarianceDecomposition')
                if isfield(oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean,name)
                    % Nothing to do.
                    return
                end
            end
        end
    end
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
            posterior_moments(tmp,1,mh_conf_sig);
    end
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean.' name ' = post_mean;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.median.' name ' = post_median;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.variance.' name ' = post_var;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdinf.' name ' = hpd_interval(1);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdsup.' name ' = hpd_interval(2);']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.deciles.' name ' = post_deciles;']);
    eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.density.' name ' = density;']);