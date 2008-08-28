function independent_metropolis_hastings(TargetFun,ProposalFun,xparam1,vv,mh_bounds,varargin)
% Independent Metropolis-Hastings algorithm.
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o varargin              list of argument following mh_bounds
%  
% OUTPUTS 
%   None
%
% ALGORITHM
%   Metropolis-Hastings.
%
% SPECIAL REQUIREMENTS
%   None.
%  
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.
global M_ options_ bayestopt_
%%%%
%%%% Initialization of the independent metropolis-hastings chains.
%%%%
[ ix2, ilogpo2, ModelName, MhDirectoryName, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns, d ] = ...
    metropolis_hastings_initialization(TargetFun,xparam1,vv,mh_bounds,varargin{:});
xparam1 = transpose(xparam1);
OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
    ProposalDensity = 'multivariate_normal_pdf';
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
    ProposalDensity = 'multivariate_student_pdf';
end
load([MhDirectoryName '/' ModelName '_mh_history'],'record');
%%%%
%%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%%%%
InitSizeArray = min([MAX_nruns*ones(nblck) nruns],[],2);
jscale = diag(bayestopt_.jscale);
for b = fblck:nblck
    if (options_.load_mh_file~=0)  & (fline(b)>1) & OpenOldFile(b)
        load(['./' MhDirectoryName '/' ModelName '_mh' int2str(NewFile(b)) ...
              '_blck' int2str(b) '.mat'])
        x2 = [x2;zeros(InitSizeArray(b)-fline(b)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(b)-fline(b)+1,1)];
        OpenOldFile(b) = 0;
    else
        x2 = zeros(InitSizeArray(b),npar);
        logpo2 = zeros(InitSizeArray(b),1);
    end
    hh = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(nblck) ')...']);
    set(hh,'Name','Metropolis-Hastings');
    isux = 0;
    jsux = 0;
    irun = fline(b);
    j = 1;
    while j <= nruns(b)
        par = feval(ProposalFun, xparam1, d * jscale, n); 
        if all(par(:)>mh_bounds(:,1)) && all(par(:)<mh_bounds(:,2))
            logpost = - feval(TargetFun,par(:),varargin{:});
        else
            logpost = -inf;
        end
        r = logpost - ilogpo2(b) + ...
            log(feval(ProposalDensity, ix2(b,:), xparam1, d, n)) - ...
            log(feval(ProposalDensity, par, xparam1, d, n));
        if (logpost > -inf) && (log(rand) < r)
            x2(irun,:) = par;
            ix2(b,:) = par;
            logpo2(irun) = logpost;
            ilogpo2(b) = logpost;
            isux = isux + 1;
            jsux = jsux + 1;
        else
            x2(irun,:) = ix2(b,:);
            logpo2(irun) = ilogpo2(b);
        end
        prtfrc = j/nruns(b);
        waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
        if (irun == InitSizeArray(b)) | (j == nruns(b)) % Now I save the simulations
            save([MhDirectoryName '/' ModelName '_mh' int2str(NewFile(b)) '_blck' int2str(b)],'x2','logpo2');
            InitSizeArray(b) = min(nruns(b)-j,MAX_nruns);
            fidlog = fopen([MhDirectoryName '/metropolis.log'],'a');
            fprintf(fidlog,['\n']);
            fprintf(fidlog,['%% Mh' int2str(NewFile(b)) 'Blck' int2str(b) ' (' datestr(now,0) ')\n']);
            fprintf(fidlog,' \n');
            fprintf(fidlog,['  Number of simulations.: ' int2str(length(logpo2)) '\n']);
            fprintf(fidlog,['  Acceptation rate......: ' num2str(jsux/length(logpo2)) '\n']);
            fprintf(fidlog,['  Posterior mean........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(mean(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(mean(logpo2)) '\n']);
            fprintf(fidlog,['  Minimum value.........:\n']);;
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(min(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(min(logpo2)) '\n']);
            fprintf(fidlog,['  Maximum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(max(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(max(logpo2)) '\n']);
            fprintf(fidlog,' \n');
            fclose(fidlog);
            jsux = 0;
            if j == nruns(b) % I record the last draw...
                record.LastParameters(b,:) = x2(end,:);
                record.LastLogLiK(b) = logpo2(end);
            end
            if InitSizeArray(b)
                x2 = zeros(InitSizeArray(b),npar);
                logpo2 = zeros(InitSizeArray(b),1);
                NewFile(b) = NewFile(b) + 1;
                irun = 0;
            else% InitSizeArray is equal to zero because we are at the end of an mc chain.
                InitSizeArray(b) = min(nruns(b),MAX_nruns);
            end
        end
        j=j+1;
        irun = irun + 1;
    end% End of the simulations for one mh-block.
    record.AcceptationRates(b) = isux/j;
    close(hh);
end% End of the loop over the mh-blocks.
record.Seeds.Normal = randn('state');
record.Seeds.Unifor = rand('state');
save([MhDirectoryName '/' ModelName '_mh_history'],'record');
disp(['MH: Number of mh files			: ' int2str(NewFile(1)) ' per block.'])
disp(['MH: Total number of generated files	: ' int2str(NewFile(1)*nblck) '.'])
disp(['MH: Total number of iterations 		: ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp(' ')