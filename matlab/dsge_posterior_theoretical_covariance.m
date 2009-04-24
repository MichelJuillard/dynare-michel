function [nvar,vartan,CovarFileNumber] = dsge_posterior_theoretical_covariance(SampleSize,M_,options_,oo_)
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

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;    

% Get informations about the _posterior_draws files.
DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
NumberOfDrawsFiles = length(DrawsFiles);

% Number of lines in posterior data files.
MaXNumberOfCovarLines = ceil(options_.MaxNumberOfBytes/(nvar*(nvar+1)/2)/8);

if SampleSize<=MaXNumberOfCovarLines
    Covariance_matrix = zeros(SampleSize,nvar*(nvar+1)/2);
    NumberOfCovarFiles = 1;
else
    Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
    NumberOfLinesInTheLastCovarFile = mod(SampleSize,MaXNumberOfCovarLines);
    NumberOfCovarFiles = ceil(SampleSize/MaXNumberOfCovarLines);
end

NumberOfCovarLines = rows(Covariance_matrix);
CovarFileNumber = 1;

% Compute 2nd order moments and save them in *_Posterior2ndOrderMoments* files
linea = 0;
for file = 1:NumberOfDrawsFiles
    load([M_.dname '/metropolis/' DrawsFiles(file).name ],'pdraws');
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
        for i=1:nvar
            for j=i:nvar
                Covariance_matrix(linea,symmetric_matrix_index(i,j,nvar)) = tmp{1}(i,j);
            end
        end
        if linea == NumberOfCovarLines
            save([ M_.dname '/metropolis/' M_.fname '_Posterior2ndOrderMoments' int2str(CovarFileNumber) '.mat' ],'Covariance_matrix');
            CovarFileNumber = CovarFileNumber + 1;
            linea = 0;
            test = CovarFileNumber-NumberOfCovarFiles;
            if ~test% Prepare the last round...
                Covariance_matrix = zeros(NumberOfLinesInTheLastCovarFile,nvar*(nvar+1)/2);
                NumberOfCovarLines = NumberOfLinesInTheLastCovarFile;
                CovarFileNumber = CovarFileNumber - 1;
            elseif test<0
                Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
            else
                clear('Covariance_matrix');
            end
        end
    end
end

options_.ar = nar;