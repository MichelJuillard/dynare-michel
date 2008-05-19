function dsge_posterior_theoretical_covariance()
% This function estimates the posterior density of the endogenous
% variables second order moments. 
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
%  
% part of DYNARE, copyright Dynare Team (2007-2008)
% Gnu Public License.

global M_ options_ oo_

type = 'posterior';% To be defined as a input argument later...
NumberOfSimulations = 800;% To be defined in a global structure...

% Set varlist (vartan)
[ivar,vartan] = set_stationary_variables_list;

% Set various parameters & Check or create files and directories
if strcmpi(type,'posterior')
    MhDirectoryName = CheckPath('metropolis');
else
    MhDirectoryName = CheckPath('prior');
end
fname = [ MhDirectoryName '/' M_.fname];
%save([fname '_Posterior2ndOrder'],'varlist');
DrawsFiles = dir([fname '_' type '_draws*' ]);
if ~rows(DrawsFiles)
    if strcmpi(type,'posterior')
        SampleAddress = selec_posterior_draws(NumberOfSimulations,1);
    else% (samples from the prior) To be done later...
    end
    DrawsFiles = dir([fname '_' type '_draws*']);
end

% Get the number of stationary endogenous variables.
nvar = length(ivar);



nar = options_.ar;% Saves size of the auto-correlation function.
options_.ar = 0;% Set the size of the auto-correlation function.

NumberOfDrawsFiles = rows(DrawsFiles);
MaXNumberOfCovarLines = ceil(options_.MaxNumberOfBytes/(nvar*(nvar+1)/2)/8);

if NumberOfSimulations<=MaXNumberOfCovarLines
    Covariance_matrix = zeros(NumberOfSimulations,nvar*(nvar+1)/2);
    NumberOfCovarFiles = 1;
else
    Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
    NumberOfLinesInTheLastCovarFile = mod(NumberOfSimulations,MaXNumberOfCovarLines);
    NumberOfCovarFiles = ceil(NumberOfSimulations/MaXNumberOfCovarLines);
end

NumberOfCovarLines = rows(Covariance_matrix);
CovarFileNumber = 1;

% Compute 2nd order moments and save them in *_Posterior2ndOrderMoments* files
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
            for j=i:nvar
                Covariance_matrix(linea,idx(i,j,nvar)) = tmp{1}(i,j);
            end
        end
        if linea == NumberOfCovarLines
            save([fname '_Posterior2ndOrderMoments' int2str(CovarFileNumber)],'Covariance_matrix');
            CovarFileNumber = CovarFileNumber + 1;
            linea = 0;
            test = CovarFileNumber-NumberOfCovarFiles;
            if ~(CovarFileNumber-NumberOfCovarFiles)% Prepare the last round...
                Covariance_matrix = zeros(NumberOfLinesInTheLastCovarFile,nvar*(nvar+1)/2);
                NumberOfCovarLines = NumberOfLinesInTheLastCovarFile;
            elseif CovarFileNumber-NumberOfCovarFiles<0;
                Covariance_matrix = zeros(MaXNumberOfCovarLines,nvar*(nvar+1)/2);
            else
                clear('Covariance_matrix');    
            end
        end
    end
end
options_.ar = nar; clear('pdraws','tmp');

% Compute statistics and save in oo_
for i=1:nvar
    for j=i:nvar
        i1 = 1;
        tmp = zeros(NumberOfSimulations,1);
        for file = 1:NumberOfDrawsFiles
            load([fname '_Posterior2ndOrderMoments' int2str(file)]);
            i2 = i1 + rows(Covariance_matrix) - 1;
            tmp(i1:i2) = Covariance_matrix(:,idx(i,j,nvar));
            i1 = i2+1;
        end
        [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = ...
            posterior_moments(tmp,1,options_.mh_conf_sig);
        name = fieldname(i,j,vartan);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.mean.' name ' = post_mean;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.median.' name ' = post_median;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.variance.' name ' = post_var;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdinf.' name ' = hpd_interval(1);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.hpdsup.' name ' = hpd_interval(2);']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.deciles.' name ' = post_deciles;']);
        eval(['oo_.PosteriorTheoreticalMoments.dsge.covariance.density.' name ' = density;']);
    end
end



    
function k = idx(i,j,n)
    k = (i-1)*n+j-i*(i-1)/2;


function r = rows(M)
    r = size(M,1);

function name = fieldname(i,j,vlist)
    n1 = deblank(vlist(i,:));
    n2 = deblank(vlist(j,:));
    name = [n1 '.' n2];