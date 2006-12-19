function metropolis(TargetFun,xparam1,vv,mh_bounds,varargin)
% Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o gend       [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data       [double]   (T*n) matrix of data ==> varargin{2}.
%  
% OUTPUTS 
%   None.  
%
%
% ALGORITHM 
%   Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
global M_ options_ bayestopt_

bayestopt_.penalty = 1e8;

MhDirectoryName = CheckPath('metropolis');

nblck = options_.mh_nblck;
nruns = ones(nblck,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
d = chol(vv);
options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);

if options_.load_mh_file == 0
    %% Here we start a new metropolis-hastings, previous draws are not
    %% considered.
    if nblck > 1
        disp('MH: Multiple chains mode.')
    else
        disp('MH: One Chain mode.')
    end
    %% Delete old mh files...
    files = dir([ MhDirectoryName '/' M_.fname '_mh*_blck*.mat']);
    if length(files)
        delete([ MhDirectoryName '/' M_.fname '_mh*_blck*.mat']);
        disp('MH: Old _mh files succesfully erased!')
    end
    %% Initial values...
    if nblck > 1% Case 1: multiple chains
        disp('MH: Searching for initial values...')
        ix2 = zeros(nblck,npar);
        ilogpo2 = zeros(nblck,1);
        for j=1:nblck
            validate	= 0;
            init_iter	= 0;
            trial = 1;
            while validate == 0 & trial <= 10 
                candidate = options_.mh_init_scale*randn(1,npar)*d + transpose(xparam1);
                if all(candidate' > mh_bounds(:,1)) & all(candidate' < mh_bounds(:,2)) 
                    ix2(j,:) = candidate;
                    ilogpo2(j) = - feval(TargetFun,ix2(j,:)',varargin{:});
                    j = j+1;
                    validate = 1;
                end
                init_iter = init_iter + 1;
                if init_iter > 100 & validate == 0
                    disp(['MH: I couldn''t get a valid initial value in 100 trials.'])
                    disp(['MH: You should Reduce mh_init_scale...'])
                    disp(sprintf('MH: Parameter mh_init_scale is equal to %f.',options_.mh_init_scale))
                    options_.mh_init_scale = input('MH: Enter a new value...  ');
                    trial = trial+1;
                end
            end
            if trial > 10 & ~validate
                disp(['MH: I''m unable to find a starting value for block ' int2str(j)])
                return
            end
        end
        disp('MH: Initial values found!')
        disp(' ')
    else% Case 2: one chain (we start from the posterior mode)
        candidate = transpose(xparam1);
        if all(candidate' > mh_bounds(:,1)) & all(candidate' < mh_bounds(:,2)) 
            ix2 = candidate;
            ilogpo2 = - feval(TargetFun,ix2',varargin{:});
            disp('MH: Initialization at the posterior mode.')
            disp(' ')
        else
            disp('MH: Initialization failed...')
            disp('MH: The posterior mode lies outside the prior bounds.')
            return
        end
    end
    fblck = 1;
    fline = ones(nblck,1);
    NewFile = ones(nblck,1);
    %%
    %% Creation of the mh-history file:
    %%
    file = dir([MhDirectoryName '/'  M_.fname '_mh_history.mat']);
    if length(files)
        delete([ MhDirectoryName '/' M_.fname '_mh_history.mat']);
        disp('MH: Old mh_history file succesfully erased!')
    end
    AnticipatedNumberOfFiles = ceil(nruns(1)/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = nruns(1) - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.Nblck = nblck;
    record.MhDraws = zeros(1,3);
    record.MhDraws(1,1) = nruns(1);
    record.MhDraws(1,2) = AnticipatedNumberOfFiles;
    record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.AcceptationRates = zeros(1,nblck);
    record.Seeds.Normal = randn('state');
    record.Seeds.Unifor = rand('state');
    record.InitialParameters = ix2;
    record.InitialLogLiK = ilogpo2;
    record.LastParameters = zeros(nblck,npar);
    record.LastLogLiK = zeros(nblck,1);
    record.LastFileNumber = AnticipatedNumberOfFiles ;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    save([MhDirectoryName '/' M_.fname '_mh_history'],'record');  
elseif options_.load_mh_file == 1
    %% Here we consider previous mh files (previous mh did not crash).
    disp('MH: I''m loading past metropolis-hastings simulations...')
    file = dir([ MhDirectoryName '/'  M_.fname '_mh_history.mat' ]);
    files = dir([ MhDirectoryName '/' M_.fname '_mh*.mat']);
    if ~length(files)
        disp('MH:: FAILURE! there is no MH file to load here!')
        return
    end
    if ~length(file)
        disp('MH:: FAILURE! there is no MH-history file!')
        return
    else
        load([ MhDirectoryName '/'  M_.fname '_mh_history'])
    end
    past_number_of_blocks = record.Nblck;
    if past_number_of_blocks ~= nblck
        disp('MH:: The specified number of blocks doesn''t match with the previous number of blocks!')
        disp(['MH:: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' ...
              int2str(past_number_of_blocks) '.'])
        disp(['MH:: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
        nblck = past_number_of_blocks;
        options_.mh_nblck = nblck;
    end
    % I read the last line of the last mh-file for initialization 
    % of the new metropolis-hastings simulations:
    LastFileNumber = record.LastFileNumber;
    LastLineNumber = record.LastLineNumber;
    if LastLineNumber < MAX_nruns
        NewFile = ones(nblck,1)*LastFileNumber;
    else
        NewFile = ones(nblck,1)*(LastFileNumber+1);
    end
    ilogpo2 = record.LastLogLiK;
    ix2 = record.LastParameters;
    fblck = 1;
    fline = ones(nblck,1)*(LastLineNumber+1);
    NumberOfPreviousSimulations = sum(record.MhDraws(:,1),1);
    record.MhDraws = [record.MhDraws;zeros(1,3)];
    NumberOfDrawsWrittenInThePastLastFile = MAX_nruns - LastLineNumber;
    NumberOfDrawsToBeSaved = nruns(1) - NumberOfDrawsWrittenInThePastLastFile;
    AnticipatedNumberOfFiles = ceil(NumberOfDrawsToBeSaved/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = NumberOfDrawsToBeSaved - (AnticipatedNumberOfFiles-1)*MAX_nruns;  
    record.LastFileNumber = LastFileNumber + AnticipatedNumberOfFiles;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MhDraws(end,1) = nruns(1);
    record.MhDraws(end,2) = AnticipatedNumberOfFiles;
    record.MhDraws(end,3) = AnticipatedNumberOfLinesInTheLastFile;
    randn('state',record.Seeds.Normal);
    rand('state',record.Seeds.Unifor);
    save([MhDirectoryName '/' M_.fname '_mh_history'],'record');
    disp(['MH: ... It''s done. I''ve loaded ' int2str(NumberOfPreviousSimulations) ' simulations.'])
    disp(' ')
elseif options_.load_mh_file == -1
    %% The previous metropolis-hastings crashed before the end! I try to
    %% recover the saved draws...
    disp('MH: Recover mode!')
    disp(' ')
    file = dir([MhDirectoryName '/'  M_.fname '_mh_history.mat']);
    if ~length(file)
        disp('MH:: FAILURE! there is no MH-history file!')
        return
    else
        load([ MhDirectoryName '/'  M_.fname '_mh_history'])
    end
    nblck = record.Nblck;
    options_.mh_nblck = nblck;
    if size(record.MhDraws,1) == 1
        OldMh = 0;% The crashed metropolis was the first session.
    else
        OldMh = 1;% The crashed metropolis wasn't the first session.
    end
    %% Default initialization:
    if OldMh
        ilogpo2 = record.LastLogLiK;
        ix2 = record.LastParameters;
    else
        ilogpo2 = record.InitialLogLiK;
        ix2 = record.InitialParameters;
    end
    %% Set "NewFile":
    if OldMh
        LastLineNumberInThePreviousMh = record.MhDraws(end-1,3);
        LastFileNumberInThePreviousMh = sum(record.MhDraws(1:end-1,2),1);
        if LastLineNumberInThePreviousMh < MAX_nruns
            NewFile = ones(nblck,1)*LastFileNumberInThePreviousMh;
        else
            NewFile = ones(nblck,1)*(LastFileNumberInThePreviousMh+1);
        end
    else
        NewFile = ones(nblck,1);
    end
    %% Set fline (First line):
    if OldMh
        fline = ones(nblck,1)*(LastLineNumberInThePreviousMh+1);
    else
        fline = ones(nblck,1);
    end
    %% Set fblck (First block):
    fblck = 1;
    %% How many mh files should we have ?
    ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
    ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*nblck;
    %% I count the total number of saved mh files...
    AllMhFiles = dir([MhDirectoryName '/' M_.fname '_mh*_blck*.mat']);
    TotalNumberOfMhFiles = length(AllMhFiles);
    %% I count the number of saved mh files per block
    NumberOfMhFilesPerBlock = zeros(nblck,1);
    for i = 1:nblck
        BlckMhFiles = dir([ MhDirectoryName '/' M_.fname '_mh*_blck' int2str(i) '.mat']);
        NumberOfMhFilesPerBlock(i) = length(BlckMhFiles);
    end
    tmp = NumberOfMhFilesPerBlock(1);
    %% Is there a chain with less mh files than expected ? 
    CrashedBlck = 1; b = 1;
    while b <= nblck
        if  NumberOfMhFilesPerBlock(b) < ExpectedNumberOfMhFilesPerBlock
            CrashedBlck = b;% YES, chain b!
            disp(['MH: Chain ' int2str(b) ' is uncomplete!'])
            break
        else
            disp(['MH: Chain ' int2str(b) ' is complete!'])
        end
        b = b+1;
    end
    if b>nblck
        disp('MH: You probably don''t need to recover a previous crashed metropolis...')
        disp('    or Dynare is unable to recover it.')
        error('I stop. You should modify your mod file...')
    end
    %% The new metropolis-hastings should start from chain... (fblck=CrashedBlck)
    fblck = CrashedBlck;
    %% How many mh-files are saved in this block ?
    NumberOfSavedMhFilesInTheCrashedBlck = NumberOfMhFilesPerBlock(CrashedBlck);
    %% How many mh-files were saved in this block during the last session
    %% (if there was a complete session before the crash)
    if OldMh
        ante = sum(record.MhDraws(1:end-1,2),1);
        NumberOfSavedMhFilesInTheCrashedBlck = NumberOfSavedMhFilesInTheCrashedBlck - ante;
    end
    %% Is the last mh-file of the previous session full ?
    %% (if there was a complete session before the crash)
    if OldMh && ~NumberOfSavedMhFilesInTheCrashedBlck
        load([MhDirectoryName '/' M_.fname '_mh' int2str(ante) '_blck' int2str(CrashedBlck) '.mat'],'logpo2');
        if length(logpo2) == MAX_nruns
            IsTheLastFileOfThePreviousMhFull = 1;
            NumberOfCompletedMhFiles = NumberOfMhFilesPerBlock(CrashedBlck);
            reste = 0;
        else
            IsTheLastFileOfThePreviousMhFull = 0;
            NumberOfCompletedMhFiles = ante-1;
            reste = MAX_nruns-LastLineNumberInThePreviousMh; 
        end
    elseif OldMh && NumberOfSavedMhFilesInTheCrashedBlck
        IsTheLastFileOfThePreviousMhFull = 1;
        NumberOfCompletedMhFiles = NumberOfMhFilesPerBlock(CrashedBlck);
        reste = 0;
    elseif ~OldMh && NumberOfSavedMhFilesInTheCrashedBlck
        IsTheLastFileOfThePreviousMhFull = 0;
        NumberOfCompletedMhFiles = NumberOfMhFilesPerBlock(CrashedBlck);
        reste = 0;
    elseif ~OldMh && NumberOfSavedMhFilesInTheCrashedBlck    
        IsTheLastFileOfThePreviousMhFull = 0;
        NumberOfCompletedMhFiles = 0;
        reste = 0;        
    end
    %% How many runs were saved ?
    NumberOfSavedDraws = MAX_nruns*NumberOfCompletedMhFiles + reste;
    %% Here is the number of draws we still need to complete the block:
    if OldMh
        nruns(CrashedBlck) = nruns(CrashedBlck)-(NumberOfSavedDraws-sum(record.MhDraws(1:end-1,1)));
    end
    %% I've got all the needed information... I can initialize the MH:
    if OldMh
        if NumberOfSavedMhFilesInTheCrashedBlck
            load([MhDirectoryName '/' M_.fname '_mh' ...
                  int2str(NumberOfCompletedMhFiles) '_blck' int2str(CrashedBlck) '.mat']);
            fline(CrashedBlck,1) = 1;
            NewFile(CrashedBlck) = NumberOfCompletedMhFiles+1;% NumberOfSavedMhFilesInTheCrashedBlck+1;
        else
            load([MhDirectoryName '/' M_.fname '_mh' ...
                  int2str(ante) '_blck' int2str(CrashedBlck) '.mat']);
            if reste
                fline(CrashedBlck,1) = length(logpo2)+1;
                NewFile(CrashedBlck) = LastFileNumberInThePreviousMh;
            else
                fline(CrashedBlck,1) = 1;
                NewFile(CrashedBlck) = LastFileNumberInThePreviousMh+1;
            end
        end
        ilogpo2(CrashedBlck) = logpo2(end);
        ix2(CrashedBlck,:) = x2(end,:);       
    end
end% of (if options_.load_mh_file == {0,1 or -1})
   %%%%
   %%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
   %%%%
InitSizeArray = min([MAX_nruns*ones(nblck) nruns],[],2);
for b = fblck:nblck
    if (options_.load_mh_file~=0)  & (fline(b)>1) & OpenOldFile(b)
        load(['./' MhDirectoryName '/' M_.fname '_mh' int2str(NewFile(b)) ...
              '_blck' int2str(b) '.mat'])
        x2 = [x2;zeros(InitSizeArray(b)-fline(b)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(b)-fline(b)+1,1)];
        OpenOldFile(b) = 0;
    else
        x2 = zeros(InitSizeArray(b),npar);
        logpo2 = zeros(InitSizeArray(b),1);
    end
    hh = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(nblck) ')...']);
    set(hh,'Name','Metropolis-Hastings')
    isux = 0;
    irun = fline(b);
    j = 1;
    while j <= nruns(b)
        par = randn(1,npar)*d;
        par = par.*bayestopt_.jscale' + ix2(b,:);  
        if all(par'>mh_bounds(:,1)) & all(par'<mh_bounds(:,2))
            logpost = - feval(TargetFun,par',varargin{:});
        else
            logpost = -inf;
        end
        if (logpost > -inf) & (log(rand) < logpost-ilogpo2(b))
            x2(irun,:) = par; 
            ix2(b,:) = par;
            logpo2(irun) = logpost; 
            ilogpo2(b) = logpost;
            isux = isux + 1;
        else    
            x2(irun,:) = ix2(b,:);
            logpo2(irun) = ilogpo2(b);
        end	
        prtfrc = j/nruns(b);
        waitbar(prtfrc,hh,[ '(' int2str(b) '/' int2str(nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/j)]);
        if (irun == InitSizeArray(b)) | (j == nruns(b)) % Now I save the simulations
            save([MhDirectoryName '/' M_.fname '_mh' int2str(NewFile(b)) '_blck' int2str(b)],'x2','logpo2');
            InitSizeArray(b) = min(nruns(b)-j,MAX_nruns);
            if j == nruns(b) % I record the last draw...
                record.LastParameters(b,:) = x2(end,:);
                record.LastLogLiK(b) = logpo2(end);
            end      
            if InitSizeArray(b)
                x2 = zeros(InitSizeArray(b),npar);
                logpo2 = zeros(InitSizeArray(b),1);
                NewFile(b) = NewFile(b) + 1;
                irun = 0;
            else % InitSizeArray is equal to zero because we are at the end of an mc chain.
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
save([MhDirectoryName '/' M_.fname '_mh_history'],'record');
disp(['MH: Number of mh files				: ' int2str(NewFile(1)) ' per block.'])
disp(['MH: Total number of generated files	: ' int2str(NewFile(1)*nblck) '.'])
disp(['MH: Total number of iterations 		: ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp(' ')