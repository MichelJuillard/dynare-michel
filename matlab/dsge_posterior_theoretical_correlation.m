function [nvar,vartan,CorrFileNumber] = dsge_posterior_theoretical_correlation(SampleSize,nar,M_,options_,oo_)
%function [nvar,vartan,CorrFileNumber] = dsge_posterior_theoretical_correlation(SampleSize,nar,M_,options_,oo_)
% This function estimates the posterior density of the endogenous
% variables second order moments. 
% 
% INPUTS 
%   SampleSize   [integer]
%   
%
%
% OUTPUTS 
%   None.
%
% SPECIAL REQUIREMENTS
%    Other matlab routines distributed with Dynare: set_stationary_variables_list.m
%                                                   CheckPath.m
%                                                   selec_posterior_draws.m                          
%                                                   set_parameters.m
%                                                   resol.m
%                                                   th_autocovariances.m    
%                                                   posterior_moments.m

% Copyright (C) 2007-2009 Dynare Team
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

type = 'posterior';
nodecomposition = 1;
    
% Set varlist (vartan)
[ivar,vartan] = set_stationary_variables_list;
nvar = length(ivar);

% Set the size of the auto-correlation function to nar.
oldnar = options_.ar;
options_.ar = nar;    

% Get informations about the _posterior_draws files.
DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
NumberOfDrawsFiles = length(DrawsFiles);

% Number of lines in posterior data files.
MaXNumberOfCorrLines = ceil(options_.MaxNumberOfBytes/(nvar*nvar*nar)/8);

if SampleSize<=MaXNumberOfCorrLines
    Correlation_array = zeros(SampleSize,nvar,nvar,nar);
    NumberOfCorrFiles = 1;
else
    Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
    NumberOfLinesInTheLastCorrFile = mod(SampleSize,MaXNumberOfCorrLines);
    NumberOfCorrFiles = ceil(SampleSize/MaXNumberOfCorrLines);
end

NumberOfCorrLines = rows(Correlation_array);
CorrFileNumber = 1;

% Compute 2nd order moments and save them in *_PosteriorCorrelations* files
linea = 0;
for file = 1:NumberOfDrawsFiles
    load([M_.dname '/metropolis/' DrawsFiles(file).name ]);
    NumberOfDraws = rows(pdraws);
    isdrsaved = columns(pdraws)-1;
    for linee = 1:NumberOfDraws
        linea = linea+1;
        if isdrsaved
            dr = pdraws{linee,2};
        else
            set_parameters(pdraws{linee,1});
            [dr,info] = resol(oo_.steady_state,0);
        end
        tmp = th_autocovariances(dr,ivar,M_,options_,nodecomposition);
        for i=1:nar
            Correlation_array(linea,:,:,i) = tmp{i+1};
        end
        if linea == NumberOfCorrLines
            save([ M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations' int2str(CorrFileNumber) '.mat' ],'Correlation_array');
            CorrFileNumber = CorrFileNumber + 1;
            linea = 0;
            test = CorrFileNumber-NumberOfCorrFiles;
            if ~test% Prepare the last round...
                Correlation_array = zeros(NumberOfLinesInTheLastCorrFile,nvar,nvar,nar);
                NumberOfCorrLines = NumberOfLinesInTheLastCorrFile;
                CorrFileNumber = CorrFileNumber - 1;
            elseif test<0;
                Correlation_array = zeros(MaXNumberOfCorrLines,nvar,nvar,nar);
            else
                clear('Correlation_array');
            end
        end
    end
end

options_.ar = oldnar;