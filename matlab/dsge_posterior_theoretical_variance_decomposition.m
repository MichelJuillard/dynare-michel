function dsge_posterior_theoretical_variance_decomposition()
% [stephane.adjemian@ens.fr][10-08-2007]
% This function estimates the posterior distribution of the variance
% decomposition of the observed endogenous variables.
% 
% INPUTS 
%   None.
%  
% OUTPUTS 
%   None.
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%     
% part of DYNARE, copyright Dynare Team (2007).
% Gnu Public License.
global M_ options_ oo_

type = 'posterior';% To be defined as a input argument later...
NumberOfSimulations = 800;% To be defined in a global structure...

% Set varlist (vartan)
[ivar,vartan] = set_stationary_variables_list();
ivar
vartan

% Set various parameters, Check or create files and directories &
% initialize arrays.
if strcmpi(type,'posterior')
    MhDirectoryName = CheckPath('metropolis');
else
    MhDirectoryName = CheckPath('prior');
end
fname = [ MhDirectoryName '/' M_.fname];
DrawsFiles = dir([fname '_' type '_draws*' ]);
if ~rows(DrawsFiles)
    if strcmpi(type,'posterior')
        SampleAddress = selec_posterior_draws(NumberOfSimulations,1);
    else% (samples from the prior) To be done later...
    end
    DrawsFiles = dir([fname '_' type '_draws*']);
end

nar = options_.ar;% Saves size of the auto-correlation function.
options_.ar = 0;% Set the size of the auto-correlation function.

nexo = M_.exo_nbr;
nvar = length(ivar);

NumberOfDrawsFiles = rows(DrawsFiles);
NumberOfSavedElementsPerSimulation = nvar*(nexo+1);
MaXNumberOfDecompLines = ceil(options_.MaxNumberOfBytes/NumberOfSavedElementsPerSimulation/8);

if NumberOfSimulations<=MaXNumberOfDecompLines
    Decomposition_array = zeros(NumberOfSimulations,nvar*(nexo+1));
    NumberOfDecompFiles = 1;
else
    Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*(nexo+1));
    NumberOfLinesInTheLastDecompFile = mod(NumberOfSimulations,MaXNumberOfDecompLines);
    NumberOfDecompFiles = ceil(NumberOfSimulations/MaXNumberOfDecompLines);
end

NumberOfDecompLines = rows(Decomposition_array);
DecompFileNumber = 1;


% Compute total variances (covariances are not saved) and variances
% implied by each structural shock.
linea = 0;
for file = 1:NumberOfDrawsFiles
    load([MhDirectoryName '/' DrawsFiles(file).name]);
    NumberOfDraws = rows(pdraws);
    for linee = 1:NumberOfDraws
        linea = linea+1;
        draw = pdraws(linee,:);
        set_parameters(draw);
        [dr,info] = resol(oo_.steady_state,0);
        tmp = th_autocovariances(dr,ivar);
        for i=1:nvar
            Decomposition_array(linea,i) = tmp{1}(i,i);
        end
        for i=1:nvar
            for j=1:nexo
                Decomposition_array(linea,nvar+(i-1)*nexo+j) = tmp{2}(i,j);
            end
        end
        if linea == NumberOfDecompLines
            save([fname '_PosteriorVarianceDecomposition' int2str(DecompFileNumber)],'Decomposition_array');
            DecompFileNumber = DecompFileNumber + 1;
            linea = 0;
            test = DecompFileNumber-NumberOfDecompFiles;
            if ~(DecompFileNumber-NumberOfDecompFiles)% Prepare the last round...
                Decomposition_array = zeros(NumberOfLinesInTheLastDecompFile,nvar*(nexo+1));
                NumberOfDecompLines = NumberOfLinesInTheLastDecompFile;
            elseif DecompFileNumber-NumberOfDecompFiles<0;
                Decomposition_array = zeros(MaXNumberOfDecompLines,nvar*(nexo+1));
            else
                clear('Decomposition_array');
            end
        end
    end
end
options_.ar = nar; clear('pdraws','tmp');

% Compute statistics and save in oo_

for i=1:nvar
    for j=1:nexo
        i1 = 1;
        tmp = zeros(NumberOfSimulations,1);
        for file = 1:NumberOfDrawsFiles
            load([fname '_PosteriorVarianceDecomposition' int2str(file)]);
            i2 = i1 + rows(Decomposition_array) - 1;
            tmp(i1:i2) = Decomposition_array(:,nvar+(i-1)*nexo+j);
            i1 = i2+1;
        end
        name = [ deblank(vartan(i,:)) '.' deblank(M_.exo_names(j,:)) ];
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
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(tmp,1);
        end
        
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.mean.' name ' = post_mean;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.median.' name ' = post_median;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.variance.' name ' = post_var;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdinf.' name ' = hpd_interval(1);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.hpdsup.' name ' = hpd_interval(2);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.deciles.' name ' = post_deciles;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.density.' name ' = density;']);
    end
end