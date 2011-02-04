function myoutput = random_walk_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% This function contain the most computationally intensive portion of code in
% random_walk_metropolis_hastings (the 'for xxx = fblck:nblck' loop). The branches in 'for'
% cycle and are completely independent than suitable to be executed in parallel way.
% 
% INPUTS 
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from random_walk_metropolis_hastings.m
%                                    function.
%   o fblck and nblck    [integer]   The Metropolis-Hastings chains.
%   o whoiam             [integer]   In concurrent programming a modality to refer to the differents thread running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' matlab, the slave matlab worker, local matlab, remote matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
% OUTPUTS
%   o myoutput  [struc]
%               If executed without parallel is the original output of 'for b =
%               fblck:nblck' otherwise a portion of it computed on a specific core or
%               remote machine. In this case:
%                               record;
%                               irun;
%                               NewFile;
%                               OutputFileName
%
% ALGORITHM 
%   Portion of Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS.
%   None.

% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code sutable to be executed in parallel on multi core or cluster machine,
% is removed from this function and placed in random_walk_metropolis_hastings_core.m funtion.
% Then the DYNARE parallel package contain a set of pairs matlab functios that can be executed in
% parallel and called name_function.m and name_function_core.m.
% In addition in the parallel package we have second set of functions used
% to manage the parallel computation.
%
% This function was the first function to be parallelized, later other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management funtions.


% Copyright (C) 2006-2010 Dynare Team
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

if nargin<4,
    whoiam=0;
end


global bayestopt_ estim_params_ options_  M_ oo_

% reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

TargetFun=myinputs.TargetFun;
ProposalFun=myinputs.ProposalFun;
xparam1=myinputs.xparam1;
vv=myinputs.vv;
mh_bounds=myinputs.mh_bounds;
ix2=myinputs.ix2;
ilogpo2=myinputs.ilogpo2;
ModelName=myinputs.ModelName;
fline=myinputs.fline;
npar=myinputs.npar;
nruns=myinputs.nruns;
NewFile=myinputs.NewFile;
MAX_nruns=myinputs.MAX_nruns;
d=myinputs.d;
InitSizeArray=myinputs.InitSizeArray;                    
record=myinputs.record;
varargin=myinputs.varargin;

% Necessary only for remote computing!
if whoiam
 Parallel=myinputs.Parallel;
 MasterName=myinputs.MasterName;
 DyMo=myinputs.DyMo;
 % initialize persistent variables in priordens()
 priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, ...
     bayestopt_.p3,bayestopt_.p4,1);
end

% (re)Set the penalty
bayestopt_.penalty = Inf;

MhDirectoryName = CheckPath('metropolis');

options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
end
% load([MhDirectoryName '/' ModelName '_mh_history.mat'],'record');
%%%%
%%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%%%%
jscale = diag(bayestopt_.jscale);

jloop=0;
if options_.use_dll==1
  
%%%TEST%%%
oldoptions_console_mode=options_.console_mode;
%options_.console_mode=1;

  if exist('OCTAVE_VERSION')
  oldoptions_console_mode=options_.console_mode;
    options_.console_mode=1;
  end
  for b = fblck:nblck,
    record.Seeds(b).Normal = randn('state');
    record.Seeds(b).Unifor = rand('state');
  end
  % calculate draws and get last line run in the last MH block sub-array
  irun = logMHMCMCposterior( xparam1, varargin{2}, mexext, fblck, nblck, nruns, vv)
  if irun<0
    error('Error in logMHMCMCposterior');
  end
  for b = fblck:nblck,
    record.Seeds(b).Normal = randn('state');
    record.Seeds(b).Unifor = rand('state');
    OutputFileName(b,:) = {[MhDirectoryName,filesep], [ModelName '_mh*_blck' int2str(b) '.mat']};
  end

  if exist('OCTAVE_VERSION')
    options_.console_mode=oldoptions_console_mode;
    options_.console_mode=1;
  end
  record.AcceptationRates=record_AcceptationRates;
  record.LastLogLiK=record_LastLogLiK;
  record.LastParameters=record_LastParameters;
  options_.console_mode=oldoptions_console_mode;
else
  for b = fblck:nblck,
    jloop=jloop+1;
    randn('state',record.Seeds(b).Normal);
    rand('state',record.Seeds(b).Unifor);
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
    if exist('OCTAVE_VERSION') || options_.console_mode
        diary off
        disp(' ')
    elseif whoiam
        %       keyboard;
        waitbarString = ['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...'];
        %       waitbarTitle=['Metropolis-Hastings ',options_.parallel(ThisMatlab).PcName];
        if options_.parallel(ThisMatlab).Local,
            waitbarTitle=['Local '];
        else
            waitbarTitle=[options_.parallel(ThisMatlab).PcName];
        end        
        fMessageStatus(0,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab), MasterName, DyMo);
    else,
        hh = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
        set(hh,'Name','Metropolis-Hastings');
        
    end
    isux = 0;
    jsux = 0;
    irun = fline(b);
    j = 1;
    load urand_1_1.csv
    load paramdraws_1_1.csv
    while j <= nruns(b)
        par = feval(ProposalFun, ix2(b,:), d * jscale, n);
        par=paramdraws_1_1(j,:);
        if all( par(:) > mh_bounds(:,1) ) & all( par(:) < mh_bounds(:,2) )
            try
                logpost = - feval(TargetFun, par(:),varargin{:});               
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        lurand=log(urand_1_1(j));
%        if (logpost > -inf) && (log(rand) < logpost-ilogpo2(b))
        if (logpost > -inf) && (lurand < logpost-ilogpo2(b))
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
        if exist('OCTAVE_VERSION') || options_.console_mode
            if mod(j, 10) == 0
                if exist('OCTAVE_VERSION')
                    printf('MH: Computing Metropolis-Hastings (chain %d/%d): %3.f%% done, acception rate: %3.f%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
                else
                    fprintf('   MH: Computing Metropolis-Hastings (chain %d/%d): %3.f \b%% done, acceptance rate: %3.f \b%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
                end
            end
            if mod(j,50)==0 & whoiam  
                %             keyboard;
                waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) '), ' sprintf('accept. %3.f%%%%', 100 * isux/j)];
                fMessageStatus(prtfrc,whoiam,waitbarString, '', options_.parallel(ThisMatlab), MasterName, DyMo)
            end
        else
            if mod(j, 3)==0 & ~whoiam
                waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
            elseif mod(j,50)==0 & whoiam,  
                %             keyboard;
                waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)];
                fMessageStatus(prtfrc,whoiam,waitbarString, waitbarTitle, options_.parallel(ThisMatlab), MasterName, DyMo)
            end
        end
        
        if (irun == InitSizeArray(b)) | (j == nruns(b)) % Now I save the simulations
            save([MhDirectoryName '/' ModelName '_mh' int2str(NewFile(b)) '_blck' int2str(b) '.mat'],'x2','logpo2');
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
            % size of next file in chain b
            InitSizeArray(b) = min(nruns(b)-j,MAX_nruns);
            % initialization of next file if necessary
            if InitSizeArray(b)
                x2 = zeros(InitSizeArray(b),npar);
                logpo2 = zeros(InitSizeArray(b),1);
                NewFile(b) = NewFile(b) + 1;
                irun = 0;
            end
        end
        j=j+1;
        irun = irun + 1;
    end% End of the simulations for one mh-block.
    record.AcceptationRates(b) = isux/j;
    if exist('OCTAVE_VERSION') || options_.console_mode
        if exist('OCTAVE_VERSION')
            printf('\n');
        else
            fprintf('\n');
        end
        diary on;
    elseif ~whoiam
        close(hh);
    end
    record.Seeds(b).Normal = randn('state');
    record.Seeds(b).Unifor = rand('state');
    OutputFileName(jloop,:) = {[MhDirectoryName,filesep], [ModelName '_mh*_blck' int2str(b) '.mat']};
  end% End of the loop over the mh-blocks.
end % if use_dll

myoutput.record = record;
myoutput.irun = irun;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;
