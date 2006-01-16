function metropolis(xparam1,vv,gend,data,rawdata,mh_bounds)
% stephane.adjemian@ens.fr [09-02-2005]
global M_ oo_ options_ bayestopt_ estim_params_

bayestopt_.penalty = 1e8;

DirectoryName = CheckPath('metropolis');

nblck = options_.mh_nblck;
nruns = ones(nblck,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
d = chol(vv);
options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);

if options_.load_mh_file == 0
  % Here we start a new metropolis-hastings, previous draws are not
  % considered.
  if nblck > 1
    disp('MH: Multiple chains mode.')
  else
    disp('MH: One Chain mode.')
  end
  % Delete old mh files...
  files = eval(['dir(''' DirectoryName '/'  M_.fname '_mh*_blck*.mat'');']);
  if length(files)
    delete([ DirectoryName '/' M_.fname '_mh*_blck*.mat']);
    disp('MH: Old _mh files succesfully erased!')
  end
  % Initial values...
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
	  ilogpo2(j) = -DsgeLikelihood(ix2(j,:)',gend,data);
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
      ilogpo2 = -DsgeLikelihood(ix2',gend,data);
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
  % Creation of the mh-history file:
  file = eval(['dir(''' DirectoryName '/'  M_.fname '_mh_history.mat'');']);
  if length(files)
    delete([ DirectoryName '/' M_.fname '_mh_history.mat']);
    disp('MH: Old mh_history file succesfully erased!')
  end
  AnticipatedNumberOfFiles = floor(nruns(1)/MAX_nruns);
  AnticipatedNumberOfLinesInTheLastFile = nruns(1) - AnticipatedNumberOfFiles*MAX_nruns;
  record.Nblck = nblck;
  record.MhDraws = zeros(1,3);
  record.MhDraws(1,1) = nruns(1);
  record.MhDraws(1,2) = AnticipatedNumberOfFiles+1;
  record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
  record.AcceptationRates = zeros(1,nblck);
  record.Seeds.Normal = randn('state');
  record.Seeds.Unifor = rand('state');
  record.InitialParameters = ix2;
  record.InitialLogLiK = ilogpo2;
  record.LastParameters = zeros(nblck,npar);
  record.LastLogLiK = zeros(nblck,1);
  record.LastFileNumber = AnticipatedNumberOfFiles+1;
  record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
  save([DirectoryName '/' M_.fname '_mh_history'],'record');  
elseif options_.load_mh_file == 1% Here we consider previous mh files (previous mh did not crash).
  disp('MH: I''m loading past metropolis-hastings simulations...')
  file = eval(['dir(''' DirectoryName '/'  M_.fname '_mh_history.mat'');']);
  files = eval(['dir(''' DirectoryName '/' M_.fname '_mh*.mat'');']);
  if ~length(files)
    disp('MH:: FAILURE! there is no MH file to load here!')
    return  
  end
  if ~length(file)
    disp('MH:: FAILURE! there is no MH-history file!')
    return
  else
    load([ DirectoryName '/'  M_.fname '_mh_history'])
  end
  past_number_of_blocks = record.Nblck;
  if past_number_of_blocks ~= nblck
    disp('MH:: The specified number of blocks doesn''t match with the previous number of blocks!')
    disp(['MH:: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' int2str(past_number_of_blocks) '.'])
    disp(['MH:: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
    nblck = past_number_of_blocks;
    options_.mh_nblck = nblck;
  end
  % I get the last line of the last mh-file for initialization 
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
  AnticipatedNumberOfFiles = floor(NumberOfDrawsToBeSaved/MAX_nruns);
  AnticipatedNumberOfLinesInTheLastFile = NumberOfDrawsToBeSaved - AnticipatedNumberOfFiles*MAX_nruns;  
  record.LastFileNumber = LastFileNumber + AnticipatedNumberOfFiles + 1;
  record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
  record.MhDraws(end,1) = nruns(1);
  record.MhDraws(end,2) = AnticipatedNumberOfFiles+1;
  record.MhDraws(end,3) = AnticipatedNumberOfLinesInTheLastFile;
  randn('state',record.Seeds.Normal);
  rand('state',record.Seeds.Unifor);
  save([DirectoryName '/' M_.fname '_mh_history'],'record');
  disp(['MH: ... It''s done. I''ve loaded ' int2str(NumberOfPreviousSimulations) ' simulations.'])
  disp(' ')
elseif options_.load_mh_file == -1% The previous metropolis-hastings
                                  % crashed before the end! I try to
                                  % recover the saved draws...
  disp('MH: Recover mode!')
  disp(' ')
  file = eval(['dir(''' DirectoryName '/'  M_.fname '_mh_history.mat'');']);
  if ~length(file)
    disp('MH:: FAILURE! there is no MH-history file!')
    return
  else
    load([ DirectoryName '/'  M_.fname '_mh_history'])
  end 
  nblck = record.Nblck;
  options_.mh_nblck = nblck;
  if size(record.MhDraws,1) == 1
    OldMh = 0;% The crashed metropolis was the first session.
  else
    OldMh = 1;% The crashed metropolis wasn't the first session.
  end
  % Default initialization:
  if OldMh
    ilogpo2 = record.LastLogLiK;
    ix2 = record.LastParameters;
  else
    ilogpo2 = record.InitialLogLiK;
    ix2 = record.InitialParameters;
  end
  % Set "NewFile":
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
  % Set fline (First line):
  if OldMh
    fline = ones(nblck,1)*(record.MhDraws(end-1,3)+1);
  else
    fline = ones(nblck,1);
  end
  % Set fblck (First block):
  fblck = 1;
  % How many mh files should we have ?
  ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
  ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*nblck;
  % I count the total number of saved mh files...
  AllMhFiles = eval(['dir(''' DirectoryName '/' M_.fname '_mh*_blck*.mat'');']);
  TotalNumberOfMhFiles = length(AllMhFiles);
  % I count the number of saved mh files per block
  NumberOfMhFilesPerBlock = zeros(nblck,1);
  for i = 1:nblck
    BlckMhFiles = eval(['dir(''' DirectoryName '/' M_.fname '_mh*_blck' int2str(i) '.mat'');']);
    NumberOfMhFilesPerBlock(i) = length(BlckMhFiles);
  end
  tmp = NumberOfMhFilesPerBlock(1); b = 1;
  % Is there a chain with less mh files than the first chain ? 
  CrashedBlck = 1; 
  while b <= nblck
    if  NumberOfMhFilesPerBlock(b) < ExpectedNumberOfMhFilesPerBlock
      CrashedBlck = b;% YES!
      break
    end
    b = b+1;
  end
  % The new metropolis-hastings should start from (fblck=CrashedBlck)
  fblck = CrashedBlck;
  % How many mh-files are saved in this block ?
  NumberOfSavedMhFilesInTheCrashedBlck = ...
      NumberOfMhFilesPerBlock(CrashedBlck);
  % How many mh-files were saved in this block during the last session
  % (if there was a complete session before the crash) ?
  if OldMh
    ante = sum(record.MhDraws(1:end-1,2),1);
    load(['./' DirectoryName '/' M_.fname '_mh' int2str(ante) '_blck' ...
	  int2str(CrashedBlck) '.mat'],'logpo2');
    if length(logpo2) == MAX_nruns
      IsTheLastFileOfThePreviousMhFull = 1;
    else
      IsTheLastFileOfThePreviousMhFull = 0;
    end
  else
    ante = 0;% Because the crashed session is the first one
    IsTheLastFileOfThePreviousMhFull = 1;
  end
  if ~IsTheLastFileOfThePreviousMhFull
    MhFileExist  = 1;
    MhFileNumber = ante;
    while MhFileExist
      MhFileNumber = MhFileNumber + 1;
      if ~exist(['./' DirectoryName '/' M_.fname '_mh' int2str(MhFileNumber) '_blck' int2str(CrashedBlck) '.mat'])
	MhFileExist = 0;
      end
    end
    % if MhFileNumber > ExpectedNumberOfMhFilesPerBlock % Peut être à déplacer plus haut...
    %  disp('MH : You are using the recover mode but the previous session');
    %  disp('     of the metropolis-hastings didn''t crash!')
    %  disp('MH : I stop and you shoud modify your mod file.')
    %  return
    % else
    NumberOfCompletedMhFiles = (MhFileNumber-1)-ante;
    % How many runs were saved ?
    if OldMh
      reste = MAX_nruns-record.MhDraws(end-1,3);
    else
      reste = 0
    end
    NumberOfSavedDraws = MAX_nruns*(NumberOfCompletedMhFiles) + reste;
    % Here is the number of draws we still need to complete the block:
    nruns(CrashedBlck) = nruns(CrashedBlck)-NumberOfSavedDraws; 
    % I initialize with the last saved mh file of the inccomplete
    % block:
    load(['./' DirectoryName '/' M_.fname '_mh' int2str(MhFileNumber-1) '_blck' int2str(CrashedBlck) '.mat']);
    ilogpo2(CrashedBlck) = logpo2(end);
    ix2(CrashedBlck,:) = x2(end,:);
    NewFile(CrashedBlck) = MhFileNumber;
    fline(CrashedBlck,1) = 1;
  else
    % creuser ce qui se passe dans ce cas !
    OpenOldFile(CrashedBlck) = 0;
    disp('Ok')
  end
end% of (if options_.load_mh_file == {0,1 or -1})
%%%%
%%%% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%%%%
InitSizeArray = min([MAX_nruns*ones(nblck) nruns],[],2);
for b = fblck:nblck
  if (options_.load_mh_file~=0)  & (fline(b)>1) & OpenOldFile(b)
    load(['./' DirectoryName '/' M_.fname '_mh' int2str(NewFile(b)) ...
	  '_blck' int2str(b) '.'])
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
      logpost = -DsgeLikelihood(par',gend,data);
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
      save([DirectoryName '/' M_.fname '_mh' int2str(NewFile(b)) '_blck' int2str(b)],'x2','logpo2');
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
save([DirectoryName '/' M_.fname '_mh_history'],'record');
disp(['MH: Number of mh files				: ' int2str(NewFile(1)) ' per block.'])
disp(['MH: Total number of generated files	: ' int2str(NewFile(1)*nblck) '.'])
disp(['MH: Total number of iterations 		: ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp(' ')
