function myoutput = random_walk_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)

% Copyright (C) 2006-2008 Dynare Team
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

struct2local(myinputs);


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
    if exist('OCTAVE_VERSION')
      diary off;
    elseif ~whoiam
      hh = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
      set(hh,'Name','Metropolis-Hastings');
    elseif whoiam
%       keyboard;
      waitbarString = ['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...'];
      waitbarTitle=['Metropolis-Hastings ',options_.parallel(ThisMatlab).PcName];
      fMessageStatus(0,b,waitbarString, waitbarTitle, options_.parallel(ThisMatlab), MasterName, DyMo);
      
    end
    isux = 0;
    jsux = 0;
    irun = fline(b);
    j = 1;
    while j <= nruns(b)
    par = feval(ProposalFun, ix2(b,:), d * jscale, n);
        if all( par(:) > mh_bounds(:,1) ) & all( par(:) < mh_bounds(:,2) )
                logpost = - feval(TargetFun, par(:),varargin{:});               
        else
            logpost = -inf;
        end
        if (logpost > -inf) & (log(rand) < logpost-ilogpo2(b))
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
        if exist('OCTAVE_VERSION')
          if mod(j, 10) == 0 & ~whoiam
            printf('MH: Computing Metropolis-Hastings (chain %d/%d): %3.f%% done, acception rate: %3.f%%\r', b, nblck, 100 * prtfrc, 100 * isux / j);
          end
        else
          if mod(j, 3)==0 & ~whoiam
            waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
          elseif mod(j,50)==0 & whoiam,  
%             keyboard;
            waitbarString = [ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)];
            fMessageStatus(prtfrc,b,waitbarString, waitbarTitle, options_.parallel(ThisMatlab), MasterName, DyMo)
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
    if exist('OCTAVE_VERSION')
      printf('\n');
      diary on;
    elseif ~whoiam
      close(hh);
    end
    record.Seeds(b).Normal = randn('state');
    record.Seeds(b).Unifor = rand('state');
    OutputFileName(jloop,:) = {[MhDirectoryName,'/'], [ModelName '_mh*_blck' int2str(b) '.mat']};
end% End of the loop over the mh-blocks.


myoutput.record = record;
myoutput.irun = irun;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;


