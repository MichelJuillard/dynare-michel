function [nvar,vartan,NumberOfDecompFiles] = ...
        dsge_posterior_theoretical_variance_decomposition(SampleSize,M_,options_,oo_)
% This function estimates the posterior distribution of the variance
% decomposition of the observed endogenous variables.
% 
% INPUTS 
%   None.
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
%     
% part of DYNARE, copyright Dynare Team (2007-2008).
% Gnu Public License.

type = 'posterior';% To be defined as a input argument later...

% Set varlist (vartan)
[ivar,vartan] = set_stationary_variables_list;
nvar = length(ivar);

% Set the size of the auto-correlation function to zero.
nar = options_.ar;
options_.ar = 0;    

% Get informations about the _posterior_draws files.
DrawsFiles = dir([M_.dname '/metropolis/' M_.fname '_' type '_draws*' ]);
NumberOfDrawsFiles = length(DrawsFiles);

nexo = M_.exo_nbr;


NumberOfDrawsFiles = rows(DrawsFiles);
NumberOfSavedElementsPerSimulation = nvar*(nexo+1);
MaXNumberOfDecompLines = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation/8);

if SampleSize<=MaXNumberOfDecompLines
    Decomposition_array = zeros(SampleSize,nvar*nexo);
    NumberOfDecompFiles = 1;
else
    Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*nexo);
    NumberOfLinesInTheLastDecompFile = mod(SampleSize,MaXNumberOfDecompLines);
    NumberOfDecompFiles = ceil(SampleSize/MaXNumberOfDecompLines);
end

NumberOfDecompLines = rows(Decomposition_array);
DecompFileNumber = 1;


% Compute total variances (covariances are not saved) and variances
% implied by each structural shock.
linea = 0;
for file = 1:NumberOfDrawsFiles
    load([M_.dname '/metropolis/' DrawsFiles(file).name ]);
    isdrsaved = cols(pdraws)-1;
    NumberOfDraws = rows(pdraws);
    for linee = 1:NumberOfDraws
        linea = linea+1;
        if isdrsaved
            dr = pdraws{linee,2};
        else
            set_parameters(pdraws{linee,1});
            [dr,info] = resol(oo_.steady_state,0);
        end
        tmp = th_autocovariances(dr,ivar,M_,options_);
        for i=1:nvar
            for j=1:nexo
                Decomposition_array(linea,(i-1)*nexo+j) = tmp{2}(i,j);
            end
        end
        if linea == NumberOfDecompLines
            save([M_.dname '/metropolis/' M_.fname '_PosteriorVarianceDecomposition' int2str(DecompFileNumber) '.mat' ],'Decomposition_array');
            DecompFileNumber = DecompFileNumber + 1;
            linea = 0;
            test = DecompFileNumber-NumberOfDecompFiles;
            if ~test% Prepare the last round...
                Decomposition_array = zeros(NumberOfLinesInTheLastDecompFile,nvar*nexo);
                NumberOfDecompLines = NumberOfLinesInTheLastDecompFile;
                DecompFileNumber = DecompFileNumber - 1;
            elseif test<0;
                Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*nexo);
            else
                clear('Decomposition_array');
            end
        end
    end
end

options_.ar = nar;