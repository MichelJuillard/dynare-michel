function metropolis(xparam1,vv,gend,data,rawdata,mh_bounds)
% stephane.adjemian@cepremap.cnrs.fr [07-31-2004]
% Adapted from an older version of metropolis.m 
global M_ oo_ options_ bayestopt_ estim_params_

TeX   	= options_.TeX;
nruns 	= options_.mh_replic;
truns 	= options_.mh_replic*options_.mh_nblck;
nblck 	= options_.mh_nblck;
nvx   	= estim_params_.nvx;
nvn   	= estim_params_.nvn;
ncx   	= estim_params_.ncx;
ncn   	= estim_params_.ncn;
np    	= estim_params_.np ;
nx    	= nvx+nvn+ncx+ncn+np;
npar  	= length(xparam1);
nvobs 	= size(options_.varobs,1);
horizon = options_.forecast;
update  = ~options_.posterior_mode_estimation;
SteP    = 1000;

%% Determine the value of MAX_nruns, MAX_nforc, MAX_nsmoo and MAX_ninno values
MaxNumberOfBytes = 1000000;%% This value should be adjusted
MAX_nruns = ceil(MaxNumberOfBytes/(npar+2)/8);
MAX_nforc = ceil(MaxNumberOfBytes/((options_.forecast+M_.maximum_lag)*length(oo_.steady_state))/8);
MAX_nsmoo = ceil(MaxNumberOfBytes/((size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic)*gend)/8);
MAX_ninno = ceil(MaxNumberOfBytes/(M_.exo_nbr*gend)/8);
MAX_nerro = ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8);
MAX_nfilt = ceil(MaxNumberOfBytes/((size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic)*gend)/8);
if options_.irf
	MAX_nirfs = ceil(MaxNumberOfBytes/(options_.irf*length(oo_.steady_state)*M_.exo_nbr)/8);
end
MAX_nthm1 = ceil(MaxNumberOfBytes/(length(oo_.steady_state)*8));
MAX_nthm2 = ceil(MaxNumberOfBytes/(length(oo_.steady_state)*length(oo_.steady_state)*8));
MAX_nthm3 = ceil(MaxNumberOfBytes/(length(oo_.steady_state)*M_.exo_nbr*8));
MAX_nthm4 = ceil(MaxNumberOfBytes/(length(oo_.steady_state)*options_.ar*8));

if options_.posterior_mode_estimation
  % If the initial condition of the metropolis-hastings algorithm
  % is around or on the posterior mode.
  d = chol(vv);
else
  % If the initial condition of the metropolis-hastings algorithm
  % is not related to the posterior mode. This matrix, used for
  % the jumping distribution, will be updated during the 
  % metropolis-hastings.
  d = eye(npar);
end

if options_.load_mh_file
  disp('MH: I''m loading past metropolis-hastings simulations...')
  files = eval(['dir(''' M_.fname '_mh*.mat'');']);
  if ~length(files)
    error('MH: FAILURE :: there is no MH file to load here!')    
  end
  bfiles = eval(['dir(''' M_.fname '_mh0_blck*.mat'');']);
  past_number_of_blocks = length(bfiles);
  lfile = size(files,1)/nblck-1;
  lfile
  number_of_simulations_per_old_file = zeros(lfile+1,1);
  if nblck == 1
    nops = 0;                 % Number Of Previous Simulations.
    MU = zeros(npar,1);       % Previous posterior mean of the parameters.
    SIGMA = zeros(npar,npar); % Previous posterior covariance matrix of the parameters. 
    tmp0 = -Inf;
    posterior_mode = zeros(1,npar);
    for file = 0:lfile
      eval(['load ' M_.fname '_mh' int2str(file)]);
      [MU,SIGMA] = recursive_moments(MU,SIGMA,x2,nops);
      number_of_simulations_per_old_file(file+1) = size(logpo2,1);
      nops = nops + number_of_simulations_per_old_file(file+1);  
      [tmp1,indx] = max(logpo2);
      if tmp1-tmp0>0
	posterior_mode = x2(indx,:);
	tmp0 = tmp1;
      end
    end
    lpost_mode = tmp0;
    if nruns
      ix2 = x2(end,:);
      ilogpo2 = logpo2(end);
    end  
    clear x2 logpo2 post2 file;
  else
    MU = zeros(npar,nblck);
    SIGMA = zeros(npar,npar,nblck);
    lpost_mode = zeros(nblck,1);
    posterior_mode = zeros(nblck,npar);
    for b = 1:nblck
      nops = 0;
      tmp0 = -Inf;
      for file = 0:lfile
	eval(['load ' M_.fname '_mh' int2str(file) '_blck' int2str(b)]);
	[MU(:,b),SIGMA(:,:,b)] = recursive_moments(MU(:,b),SIGMA(:,:,b),x2,nops);
	number_of_simulations_per_old_file(file+1) = size(logpo2,1);
	nops = nops + number_of_simulations_per_old_file(file+1);
	[tmp1,indx] = max(logpo2);
	if tmp1-tmp0>0
	  posterior_mode(b,:) = x2(indx,:);
	  tmp0 = tmp1;
	end
      end
      lpost_mode(b) = tmp0;	
      if nruns
	ix2(b,:) = x2(end,:); 	
	ilogpo2(b) = logpo2(end);
      end
    end
    clear x2 logpo2 post2 b file;
  end 
  disp(['MH: ... It''s done. I''ve loaded ' int2str(nops) 'simulations.'])
  disp(' ')  
end
%
%%
%%% METROPOLIS-HASTINGS:
%%
%
options_.lik_algo = 1;
if nruns
  if ~options_.load_mh_file
    MU = zeros(npar,nblck);
    SIGMA = repmat(zeros(npar,npar),[1 1 nblck]);
    % Delete old mh files...
    if nblck > 1
      disp('MH: Multiple chains mode.')
    else
      disp('MH: One Chain mode.')
    end
    files = eval(['dir(''' M_.fname '_mh*.mat'');']);
    if length(files)
      delete([M_.fname '_mh*.mat']);
      disp('MH: Old _mh files succesfully erased!')
    end
    nops = 0;   % Number Of Past Simulations.
    lfile = -1;	% Index for the last mh file.
    if nblck > 1
      disp('MH: Searching for initial values...')
      ix2 = zeros(nblck,npar);
      ilogpo2 = zeros(nblck,1);
      tmp = options_.mh_init_scale*ones(nblck,1);
      for b=1:nblck
	validate = 0;
	init_iter = 0;
	trial = 1;
	while validate == 0 & trial <= 999 
	  candidate = tmp(b)*randn(1,npar)*d + xparam1';
	  if all(candidate' > mh_bounds(:,1)) & all(candidate' < mh_bounds(:,2)) 
	    ix2(b,:) = candidate;
	    ilogpo2(b) = -DsgeLikelihood(ix2(b,:)',gend,data);
	    % b = b+1;
	    validate = 1;
	  end
	  init_iter = init_iter + 1;
	  if init_iter > 100 & validate == 0
	    disp(['MH: I couldn''t get a valid initial value in 100 trials' ...
		  ' (block' int2str(b) ')' ])
	    disp(sprintf('MH: I reduce mh_init_scale to %f',tmp(b)/1.1))
	    tmp(b) = tmp(b)/1.1; disp(' ')
	    trial = trial+1;
	  end
	end
	if trial > 999 & ~validate
	  error(['MH: I''m unable to find a starting value for block ' int2str(b)]);
	end
      end
      disp('MH: Initial values found!')
      disp(' ')
      if options_.posterior_mode_estimation
	d = repmat(d,[1,1,nblck]);
      else
	d = zeros(npar,npar,nblck);
	for b = 1:nblck
	  d(:,:,b) = eye(npar)*tmp(b)*0.05; % chol is useless here since SIGMA is diagonal
	end  
      end
      clear tmp;
    else
      candidate = transpose(xparam1);
      if all(transpose(candidate) > mh_bounds(:,1)) & all(transpose(candidate) < mh_bounds(:,2)) 
	ix2 = candidate;
	ilogpo2 = -DsgeLikelihood(transpose(ix2),gend,data);
	if options_.posterior_mode_estimation
	  disp('MH: Initialization at the posterior mode.')
	else
	  disp('MH: Initialization at the prior mode.')
	end   
	  disp(' ')
      else
	disp('MH: Initialization failed...')
	if options_.posterior_mode_estimation
	  error('MH: The posterior mode lies outside the prior bounds.')
	end  
      end
    end
  else
    if length(bfiles)>0 & (past_number_of_blocks-nblck)
      disp('MH: The specified number of blocks doesn''t match with the previous number of blocks!')
      disp(['MH: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' ...
	    int2str(past_number_of_blocks) '.'])
      disp(['MH: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
      nblck = past_number_of_blocks;
      options_.mh_nblck = nblck;
    end
    if options_.posterior_mode_estimation
      d = repmat(d,[1,1,nblck]);
    else
      d = zeros(npar,npar,nblck);
      for b = 1:nblck
	d(:,:,b) = chol(SIGMA(:,:,b)); 
      end
    end
  end
  number_of_simulations_per_file = zeros(ceil(nruns/MAX_nruns),1);
  if (nblck-1)
    disp('Acceptation rates :')
  end
  for b = 1:nblck
    isux = 0;
    Lfile = lfile;
    hh = waitbar(0,'Please wait... Metropolis-Hastings...');
    if ~(nblck-1)
      set(hh,'Name','Metropolis-Hastings')
    else
      set(hh,'Name',['Metropolis-Hastings, Block ',int2str(b)])
    end
    if nruns <= MAX_nruns
      x2 = zeros(nruns,npar);
      x2(1,:) = ix2(b,:);
      logpo2 = zeros(nruns,1);
      logpo2(1) = ilogpo2(b);
    else
      x2 = zeros(MAX_nruns,npar);
      x2(1,:) = ix2(b,:);
      logpo2 = zeros(MAX_nruns,1);
      logpo2(1) = ilogpo2(b); 
    end
    irun = ~options_.load_mh_file; % irun=0 means that previous files have been loaded.
    rruns = nruns-irun;
    j=1;
    while j<=rruns
      irun = irun + 1;
      if irun <= MAX_nruns
	% [1] I obtain a new vector of parameters from the jumping distribution. 
	par = randn(1,npar)*d(:,:,b);
	par = par.*bayestopt_.jscale' + ix2(b,:);
	% [2] Do I keep this vector ?...
	if all( par' > mh_bounds(:,1)) & all( par' < mh_bounds(:,2))
	  logpost = -DsgeLikelihood(par',gend,data);
	else
	  logpost = -inf;
	end
	if logpost > -inf & log(rand) < logpost - ilogpo2(b) % ...Yes! The vector is saved in x2.
	  x2(irun,:) = par; 
	  ix2(b,:) = par;
	  logpo2(irun) = logpost; 
	  ilogpo2(b) = logpost;
	  isux = isux + 1;
	else % No! Vector "par" is discarded and the past realization duplicated in x2.    
	  x2(irun,:) = ix2(b,:);
	  logpo2(irun) = ilogpo2(b);
	end
	prtfrc = j/nruns;
	waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));
	% [3] I update the mean and variance of the parameters
	if update
	  OldMU(:,b) = MU(:,b);
	  OldSIGMA(:,:,b) = SIGMA(:,:,b);
	  t = nops+j;
	  MU(:,b) = OldMU(:,b) + (x2(irun,:)'-OldMU(:,b))/t;
	  SIGMA(:,:,b) = OldSIGMA(:,:,b) + (x2(irun,:)'*x2(irun,:)-MU(:,b)*MU(:,b)'-OldSIGMA(:,:,b))/t + ...
	      (t-1)*(OldMU(:,b)*OldMU(:,b)' - MU(:,b)*MU(:,b)')/t;
	  UpdateTheJumpingDistribution = ~(floor(j/SteP)<(j/SteP));
	  if UpdateTheJumpingDistribution 
	    d(:,:,b) = chol(SIGMA(:,:,b));  
	  end
	end  
      else
	post2 = exp(logpo2);
	if ~(nblck-1)
	  save([M_.fname '_mh' int2str(Lfile+1)],'x2','logpo2','post2');
	else
	  save([M_.fname '_mh' int2str(Lfile+1) '_blck' int2str(b)],'x2','logpo2','post2');
	end
	number_of_simulations_per_file(Lfile+1-lfile,1) = MAX_nruns;
	clear x2 logpo2 post2;
	x2 = zeros(MAX_nruns,npar);
	logpo2 = zeros(MAX_nruns,1);
	Lfile = Lfile + 1;
	irun = 0;
	j = j - 1;
      end
      j = j + 1;
    end
    if nruns <= MAX_nruns
      post2 = exp(logpo2);
      if ~(nblck-1)
	save([M_.fname '_mh' int2str(lfile+1)],'x2','logpo2','post2');
      else
	save([M_.fname '_mh' int2str(lfile+1) '_blck' int2str(b)],'x2','logpo2','post2');	
      end
      number_of_simulations_per_file = nruns;
      clear post2 x2 logpo2;
    elseif irun <= MAX_nruns
      x2 = x2(1:irun,:);
      logpo2 = logpo2(1:irun,1); 
      post2 = exp(logpo2);
      if ~(nblck-1)
	save([M_.fname '_mh' int2str(Lfile+1)],'x2','logpo2','post2');
      else
	save([M_.fname '_mh' int2str(Lfile+1) '_blck' int2str(b)],'x2','logpo2','post2');
	clear post2 x2 logpo2;
      end      
      number_of_simulations_per_file(Lfile+1-lfile,1) = irun;
    end
    close(hh)
    if ~(nblck-1)
      disp(sprintf('Acceptation rate : %f',isux/nruns))
    else
      disp(sprintf('Block %d: %f',b,isux/nruns))      
    end
  end
  disp(' ')
  disp(['MH: Total number of iterations : ' int2str(nops+nruns) '.'])
end
if options_.load_mh_file & nruns
  number_of_simulations_per_file = cat(1,number_of_simulations_per_old_file, ...
				       number_of_simulations_per_file);
elseif options_.load_mh_file & ~nruns
  number_of_simulations_per_file = number_of_simulations_per_old_file;
end
nfile = size(number_of_simulations_per_file,1)-1;
cumulated_number_of_simulations_per_file = cumsum(number_of_simulations_per_file);
disp(['MH: Number of mh files			: ' int2str(nfile+1) ' per block.'])
disp(['MH: Total number of generated files	: ' int2str((nfile+1)*nblck) '.'])
disp(['MH: Total number of iterations 		: ' int2str(nops+nruns) '.'])
disp('MH: Number of simulations per file: ')
for i=0:nfile
	disp(sprintf('    The number of simulations in file %d is: %d.',i,number_of_simulations_per_file(i+1)))
end
disp(' ')
nsim = nops+nruns;
%
%%
%%%
%%%%
%%%%% MCMC convergence diagnostics
%%%%
%%%
%%
%
if ~options_.nodiagnostic & nblck > 1
  McMCDiagnostics(1000,0.2,npar,nsim,lfile,nblck, ...
		  number_of_simulations_per_file,cumulated_number_of_simulations_per_file);  
end  
%%
%% Now i discard some simulations...
%%
trun = cumulated_number_of_simulations_per_file(nfile+1);
irun = floor(options_.mh_drop*trun)+1;
ffil = 0;       % The first MH file we have to read...
ifil = irun;    % and the first line we have to read in this file.
for ffil = 0:nfile
	if irun <= cumulated_number_of_simulations_per_file(ffil+1)
    	break
    end
    ifil = ifil-number_of_simulations_per_file(ffil+1);
end
trun = trun-irun+1;
fprintf('MH: I''ll use mh-files %d to %d.\n',ffil,nfile);
fprintf('MH: In mh-file number %d i''ll start at line %d.\n',ffil,ifil);
fprintf('MH: Finally the total number of simulations is %d.\n',trun);
disp(' ');
%
%%
%%%
%%%%
%%%%% Modified harmonic mean
%%%%
%%%
%%
%
fprintf('MH: I''m computing the posterior mean... ');
if nblck > 1
  MU = zeros(1,npar);
  lpost_mode = -Inf;
  for  b = 1:nblck
    instr = [M_.fname '_mh' int2str(ffil) '_blck' int2str(b)];
    eval(['load ' instr]); clear post2;
    MU(1,:) = MU(1,:) + sum(x2(ifil:end,:),1);
    lpost_mode = max(lpost_mode,max(logpo2(ifil:end,1)));
  end
  for n = ffil+1:nfile
    for b = 1:nblck
      instr = [M_.fname '_mh' int2str(n) '_blck' int2str(b)];
      eval(['load ' instr]);
      clear post2;
      MU(1,:) = MU(1,:) + sum(x2,1);
      lpost_mode = max(lpost_mode,max(logpo2));
    end
  end
  clear x2 logpo2;
else
  MU = zeros(1,npar);
  lpost_mode = -Inf;
  instr = [M_.fname '_mh' int2str(ffil)];
  eval(['load ' instr]);
  clear post2;
  MU(1,:) = MU(1,:) + sum(x2(ifil:end,:),1);
  lpost_mode = max(lpost_mode,max(logpo2(ifil:end)));
  for n = ffil+1:nfile
    instr = [M_.fname '_mh' int2str(n)];
    eval(['load ' instr]);
    clear post2;
    MU(1,:) = MU(1,:) + sum(x2,1);
    lpost_mode = max(lpost_mode,max(logpo2));
  end
  clear x2 logpo2;
end %	<=== Mean of the parameters (ok!)
MU = MU/(trun*nblck);
%	lpost_mode is the value of the log posterior kernel at the mode.	
fprintf(' Done!\n');
fprintf('MH: I''m computing the posterior covariance matrix... ');
SIGMA = zeros(npar,npar);
if nblck > 1
  for b = 1:nblck
    instr = [M_.fname '_mh' int2str(ffil) '_blck' int2str(b)];
    eval(['load ' instr]);
    clear post2 logpo2;
    SIGMA = SIGMA + transpose(x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU)*...
	    (x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU);
  end				
  for n = ffil+1:nfile
    for b = 1:nblck
      instr = [M_.fname '_mh' int2str(n) '_blck' int2str(b)];
      eval(['load ' instr]);
      clear post2 logpo2;
      SIGMA = SIGMA + transpose(x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU)*...
	      (x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU);
    end				
  end
  clear x2;
else
  instr = [M_.fname '_mh' int2str(ffil)];
  eval(['load ' instr]);
  clear post2 logpo2;
  SIGMA = SIGMA + transpose(x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU)*...
	  (x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU);
  for n = ffil+1:nfile
    instr = [M_.fname '_mh' int2str(n)];
    eval(['load ' instr]);
    clear post2 logpo2;
    SIGMA = SIGMA + transpose(x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU)*...
	    (x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU);
  end
  clear x2;
end
SIGMA =  SIGMA/(trun*nblck);%<=== Variance of the parameters (ok!)
fprintf(' Done!\n');
disp(' ');
disp('MH: I''m computing the posterior log marginale density (modified harmonic mean)... ');
detSIGMA = det(SIGMA);
invSIGMA = inv(SIGMA);
marginal = zeros(9,2);
linee = 0;
check_coverage  = 1;
increase        = 1;
while check_coverage
  for p = 0.1:0.1:0.9;
    critval = qchisq(p,npar);
    tmp = 0;
    if nblck == 1
      instr = [M_.fname '_mh' int2str(ffil)];
      eval(['load ' instr]);
      clear post2;
      EndOfFile = number_of_simulations_per_file(ffil+1);
      for i = ifil:EndOfFile
	deviation  = (x2(i,:)-MU)*invSIGMA*transpose(x2(i,:)-MU);
	if deviation <= critval
	  lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
	  tmp = tmp + exp(lftheta - logpo2(i)+lpost_mode);
	end
      end
      for k = ffil+1:nfile
	instr = [M_.fname '_mh' int2str(k)];
	eval(['load ' instr]);
	clear post2;
	EndOfFile = number_of_simulations_per_file(k+1);
	for i = 1:EndOfFile
	  deviation  = (x2(i,:)-MU)*invSIGMA*transpose(x2(i,:)-MU);
	  if deviation <= critval
	    lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
	    tmp = tmp + exp(lftheta - logpo2(i)+lpost_mode);
	  end
	end
      end
      clear x2 logpo2;
    else	
      inst = [M_.fname '_mh' int2str(ffil)];
      EndOfFile = number_of_simulations_per_file(ffil+1);
      for b=1:nblck
	instr = [inst '_blck' int2str(b)];
	eval(['load ' instr]);
	clear post2;
	for i = ifil:EndOfFile
	  for j=1:nblck
	    deviation  = (x2(i,:)-MU)*invSIGMA*transpose(x2(i,:)-MU);
	    if deviation <= critval
	      lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
	      tmp = tmp + exp(lftheta - logpo2(i)+lpost_mode);
	    end
	  end
	end
      end	
      for k = ffil+1:nfile
	inst = [M_.fname '_mh' int2str(k)];
	EndOfFile = number_of_simulations_per_file(k+1);
	for b=1:nblck
	  instr = [inst '_blck' int2str(b)];
	  eval(['load ' instr]);
	  clear post2;
	  for i = ifil:EndOfFile
	    for j=1:nblck
	      deviation  = (x2(i,:)-MU)*invSIGMA*transpose(x2(i,:)-MU);
	      if deviation <= critval
		lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
		tmp = tmp + exp(lftheta - logpo2(i)+lpost_mode);
	      end
	    end
	  end
	end	
      end
      clear x2 logpo2;
    end
    linee = linee + 1;
    warning off all
    marginal(linee,:) = [p,lpost_mode-log(tmp/(trun*nblck))];
    warning on all
  end
  if abs((marginal(9,2)-marginal(1,2))/marginal(9,2)) > 0.01 | isinf(marginal(1,2))
    if increase == 1
      disp('MH: The support of the weighting density function is not large enough...')
      disp('MH: I increase the variance of this distribution.')
      increase = 1.2*increase;
      invSIGMA = inv(SIGMA*increase);
      detSIGMA = det(SIGMA*increase);
      linee    = 0;   
    else
      disp('MH: Let me try again.')
      increase = 1.2*increase;
      invSIGMA = inv(SIGMA*increase);
      detSIGMA = det(SIGMA*increase);
      linee    = 0;
      if increase > 20
	check_coverage = 0;
	clear invSIGMA detSIGMA increase;
	disp('MH: There''s probably a problem with the modified harmonic mean estimator.')    
      end    
    end    
  else
    check_coverage = 0;
    clear invSIGMA detSIGMA increase;
    disp('MH: Modified harmonic mean estimator, done!')
  end    
end
%
%%
%%%
%%%%
%%%%% Highest Probability Intervals (coverage is given by options_.mh_conf_sig)
%%%%
%%%
%%
%
disp(' ')
fprintf('MH: I''m computing the Highest Probability Intervals... ');
post_mean = transpose(MU); clear MU;
n	= trun*nblck;
n1	= round((1-options_.mh_conf_sig)*n);
k	= zeros(n1,1);
tmp = zeros(n,1);
if nblck == 1
	for i = 1:npar
		EndOfFile = number_of_simulations_per_file(ffil+1)-ifil+1;
		instr = [M_.fname '_mh' int2str(ffil)];
		eval(['load ' instr]);
		clear post2 logpo2;
		tmp(1:EndOfFile) = x2(ifil:end,i);
		OldEndOfFile = EndOfFile;
		for f = ffil+1:nfile
	  		NewEndOfFile = number_of_simulations_per_file(f+1);
	  		instr = [M_.fname '_mh' int2str(f)];
	  		eval(['load ' instr]);
	  		clear post2 logpo2;
	  		tmp(OldEndOfFile+1:OldEndOfFile+NewEndOfFile) = x2(:,i);
	  		OldEndOfFile = OldEndOfFile + NewEndOfFile;
		end
		clear x2;
		tmp = sort(tmp);
		j2 = n-n1;
		for j1 = 1:n1
	  		k(j1) = tmp(j2)-tmp(j1);
	  		j2 = j2 + 1;
		end
		[kmin,k1] = min(k);
		min_interval(i,:) = [tmp(k1) tmp(k1)+kmin];
	end
    clear tmp;
else
	for i = 1:npar
		EndOfFile = number_of_simulations_per_file(ffil+1)-ifil+1;
		NewStartLine = 0;
		inst = [M_.fname '_mh' int2str(ffil)];
		for b = 1:nblck
	  		instr = [inst '_blck' int2str(b)];
	  		eval(['load ' instr]);
	  		clear post2 logpo2;
	  		tmp(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(ifil:end,i);
	  		NewStartLine = NewStartLine + EndOfFile;
		end
		for f = ffil+1:nfile
	  		EndOfFile = number_of_simulations_per_file(f+1);
	  		inst = [M_.fname '_mh' int2str(f)];
	  		for B = 1:nblck
	    		instr = [inst '_blck' int2str(b)];
	    		eval(['load ' instr]);
	    		clear post2 logpo2;
	    		tmp(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(:,i);
	    		NewStartLine = NewStartLine + EndOfFile;
	  		end
		end
		clear x2;
		tmp = sort(tmp);
		j2 = n-n1;
		for j1 = 1:n1
	 		k(j1) = tmp(j2)-tmp(j1);
	  		j2 = j2 + 1;
		end
		[kmin,k1] = min(k);
		min_interval(i,:) = [tmp(k1) tmp(k1)+kmin];
	end
    clear tmp;
end
fprintf(' Done!\n');
%
%%
%%%
%%%%
%%%%% Print results
%%%%
%%%
%%
%%
%% [1] On screen
%%
disp(' ');
disp(' ')
marginal
disp(' ')
disp(' ')
disp('ESTIMATION RESULTS')
disp(' ')
disp(sprintf('Log data density is %f.',mean(marginal(:,2))))
oo_.MarginalDensity.ModifiedHarmonicMean = mean(marginal(:,2));
pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
tit2 = sprintf('%10s %7s %10s %14s %4s %6s\n',' ','prior mean', ...
	'post. mean','conf. interval','prior','pstdev');
ip = nvx+nvn+ncx+ncn+1;
if np
	disp(' ')
    disp('parameters')
    disp(tit2)
    for i=1:np
		disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		     deblank(estim_params_.user_param_names(i,:)), ...
		     bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
		     pnames(bayestopt_.pshape(ip)+1,:), ...
		     bayestopt_.pstdev(ip)));
	   	eval(['oo_.posterior_mean.parameters.' deblank(estim_params_.user_param_names(i,:)) ' = post_mean(ip);']);
	   	eval(['oo_.posterior_hpdinf.parameters.' deblank(estim_params_.user_param_names(i,:)) ' = min_interval(ip,1);']); 
 		eval(['oo_.posterior_hpdsup.parameters.' deblank(estim_params_.user_param_names(i,:)) ' = min_interval(ip,2);']);
		ip = ip+1;
    end
end
if nvx
	ip = 1;
    disp(' ')
    disp('standard deviation of shocks')
    disp(tit2)
    for i=1:nvx
		k = estim_params_.var_exo(i,1);
		disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		     deblank(M_.exo_names(k,:)),bayestopt_.pmean(ip),post_mean(ip), ...
		     min_interval(ip,:),pnames(bayestopt_.pshape(ip)+1,:), ...
		     bayestopt_.pstdev(ip))); 
		M_.Sigma_e(k,k) = post_mean(ip)*post_mean(ip);
	   	eval(['oo_.posterior_mean.shocks_std.' deblank(M_.exo_names(k,:)) ' = post_mean(ip);']);
	   	eval(['oo_.posterior_hpdinf.shocks_std.' deblank(M_.exo_names(k,:)) ' = min_interval(ip,1);']); 
 		eval(['oo_.posterior_hpdsup.shocks_std.' deblank(M_.exo_names(k,:)) ' = min_interval(ip,2);']);
		ip = ip+1;
    end
end
if nvn
	disp(' ')
    disp('standard deviation of measurement errors')
    disp(tit2)
    ip = nvx+1;
    for i=1:nvn
		disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		     deblank(options_.varobs(estim_params_.var_endo(i,1),:)),...
		     bayestopt_.pmean(ip), ...
		     post_mean(ip),min_interval(ip,:), ...
		     pnames(bayestopt_.pshape(ip)+1,:), ...
		     bayestopt_.pstdev(ip)));
	   	eval(['oo_.posterior_mean.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = post_mean(ip);']);
	   	eval(['oo_.posterior_hpdinf.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = min_interval(ip,1);']); 
 		eval(['oo_.posterior_hpdsup.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = min_interval(ip,2);']);		      
		ip = ip+1;
    end
end
if ncx
	disp(' ')
	disp('correlation of shocks')
    disp(tit2)
    ip = nvx+nvn+1;
    for i=1:ncx
		k1 = estim_params_.corrx(i,1);
		k2 = estim_params_.corrx(i,2);
		name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
		disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		     bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
		     pnames(bayestopt_.pshape(ip)+1,:), ...
		     bayestopt_.pstdev(ip)));
	   	eval(['oo_.posterior_mean.shocks_corr.' deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:)) ' = post_mean(ip);']);
	   	eval(['oo_.posterior_hpdinf.shocks_corr.' deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:)) ' = min_interval(ip,1);']); 
 		eval(['oo_.posterior_hpdsup.shocks_corr.' deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:)) ' = min_interval(ip,2);']);      
		M_.Sigma_e(k1,k2) = post_mean(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
		M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
		ip = ip+1;
	end
end
if ncn
	disp(' ')
    disp('correlation of measurement errors')
    disp(tit2)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
		k1 = estim_params_.corrn(i,1);
		k2 = estim_params_.corrn(i,2);
		name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
		disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		     bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
		     pnames(bayestopt_.pshape(ip)+1,:), ...
		     bayestopt_.pstdev(ip))); 
	   	eval(['oo_.posterior_mean.measurement_errors_corr.' deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:)) ' = post_mean(ip);']);
	   	eval(['oo_.posterior_hpdinf.measurement_errors_corr.' deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:)) ' = min_interval(ip,1);']); 
 		eval(['oo_.posterior_hpdsup.measurement_errors_corr.' deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:)) ' = min_interval(ip,2);']);      
		ip = ip+1;
	end
end
%%
%% [1] In a TeX file
%%
if TeX 
	if np
		ip = nvx+nvn+ncx+ncn+1;
		fidTeX = fopen([M_.fname '_MH_Posterior_1.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (parameters)\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:np
  			fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
	  			deblank(estim_params_.tex(i,:)), ...
	  			deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
	  			bayestopt_.pmean(ip), ...
	  			bayestopt_.pstdev(ip), ...
	  			post_mean(ip), ...
	  			min_interval(ip,1), ...
	  			min_interval(ip,2));
  			ip = ip+1;
		end   
		fprintf(fidTeX,'\\hline\\hline \n');
		fprintf(fidTeX,'\\end{tabular}\n ');    
		fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (parameters)}\n ');
		fprintf(fidTeX,'\\label{Table:MhPosterior:1}\n');
		fprintf(fidTeX,'\\end{table}\n');
		fprintf(fidTeX,'} \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
  	end
	if nvx
		ip = 1;
		fidTeX = fopen([M_.fname '_MH_Posterior_2.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (standard deviation of structural shocks)\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvx
  			k = estim_params_.var_exo(i,1);
  			fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
	  			deblank(M_.exo_names_tex(k,:)),...
	  			deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
	  			bayestopt_.pmean(ip), ...
	  			bayestopt_.pstdev(ip), ...
	  			post_mean(ip), ...
	  			min_interval(ip,1), ...
	  			min_interval(ip,1));
  			ip = ip+1;
		end
		fprintf(fidTeX,'\\hline\\hline \n');
		fprintf(fidTeX,'\\end{tabular}\n ');    
		fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (standard deviation of structural shocks)}\n ');
		fprintf(fidTeX,'\\label{Table:MhPosterior:2}\n');
		fprintf(fidTeX,'\\end{table}\n');
		fprintf(fidTeX,'} \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
	end
	if nvn
		ip = nvx+1;
		fidTeX = fopen([M_.fname '_MH_Posterior_3.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (standard deviation of measurement errors)\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvn
			fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
				deblank(options_.varobs_TeX(estim_params_.var_endo(i,1),:)), ...
				deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
				bayestopt_.pmean(ip), ...
				bayestopt_.pstdev(ip), ...
				post_mean(ip), ...
				min_interval(ip,1), ...
				min_interval(ip,2));
			p = ip+1;
		end
		fprintf(fidTeX,'\\hline\\hline \n');
		fprintf(fidTeX,'\\end{tabular}\n ');    
		fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (standard deviation of measurement errors)}\n ');
		fprintf(fidTeX,'\\label{Table:MhPosterior:3}\n');
		fprintf(fidTeX,'\\end{table}\n');
		fprintf(fidTeX,'} \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
	end
  	if ncx
		ip = nvx+nvn+1;
		fidTeX = fopen([M_.fname '_MH_Posterior_4.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (correlation of structural shocks)\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:ncx
			k1 = estim_params_.corrx(i,1);
		  	k2 = estim_params_.corrx(i,2);
		  	name = [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))];
		  	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
				name, ...
				deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
				bayestopt_.pmean(ip), ...
				bayestopt_.pstdev(ip), ...
				post_mean(ip), ...
				min_interval(ip,1), ...
				min_interval(ip,2));
		  	ip = ip+1;
		end
		fprintf(fidTeX,'\\hline\\hline \n');
		fprintf(fidTeX,'\\end{tabular}\n ');    
		fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (correlation of structural shocks)}\n ');
		fprintf(fidTeX,'\\label{Table:MhPosterior:4}\n');
		fprintf(fidTeX,'\\end{table}\n');
		fprintf(fidTeX,'} \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
	end
  	if ncn
		ip = nvx+nvn+ncx+1;
		fidTeX = fopen([M_.fname '_MH_Posterior_5.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (correlation of measurement errors)\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:ncn
			k1 = estim_params_.corrn(i,1);
			k2 = estim_params_.corrn(i,2);
			name = [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))];
			fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
				name, ...
				deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
				bayestopt_.pmean(ip), ...
				bayestopt_.pstdev(ip), ...
				post_mean(ip), ...
				min_interval(ip,1), ...
				min_interval(ip,2));
			ip = ip+1;
		end
		fprintf(fidTeX,'\\hline\\hline \n');
		fprintf(fidTeX,'\\end{tabular}\n ');    
		fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (correlation of structural shocks)}\n ');
		fprintf(fidTeX,'\\label{Table:MhPosterior:5}\n');
		fprintf(fidTeX,'\\end{table}\n');
		fprintf(fidTeX,'} \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
	end
end % if TeX
%                                                               
%%                                                              
%%%                                                             
%%%%                                                            
%%%%% Plot posterior distributions
%%%%                                                            
%%%                                                             
%%                                                              
%                                                               
figurename = 'Priors and posteriors';
if TeX    
	fidTeX = fopen([M_.fname '_PriorsAndPosteriors.TeX'],'w');
	fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
	fprintf(fidTeX,' \n');
end
[nbplt,nr,nc,lr,lc,nstar] = pltorg(npar);
if nbplt == 1
	h1 = figure('Name',figurename);
	if TeX
		TeXNAMES = [];
  	end
  	NAMES = []; 
  	for i=1:npar
		[borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
    		posterior_distribution(i,nfile,ffil,ifil,...
			   nblck,n,number_of_simulations_per_file,TeX);
	   	eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
	   	eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']); 
 		if TeX
  			TeXNAMES = strvcat(TeXNAMES,texnam);
		end    
		NAMES = strvcat(NAMES,nam);
		subplot(nr,nc,i);
		hh = plot(x2,f2,'-k','linewidth',2);
		set(hh,'color',[0.7 0.7 0.7]);
		hold on;
		plot(x1,f1,'-k','linewidth',2);
		plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
		box on;
		axis([borneinf bornesup 0 1.1*top]);
		title(nam,'Interpreter','none');
		hold off;
		drawnow
  	end
	eval(['print -depsc2 ' M_.fname '_PriorsAndPosteriors' int2str(1)]);
	eval(['print -dpdf ' M_.fname '_PriorsAndPosteriors' int2str(1)]);
	saveas(h1,[M_.fname '_PriorsAndPosteriors' int2str(1) '.fig']);
	if TeX
		fprintf(fidTeX,'\\begin{figure}[H]\n');
		for jj = 1:npar
			fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
		end    
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',M_.fname,int2str(1));
		fprintf(fidTeX,'\\caption{Priors and posteriors.}');
		fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(1));
		fprintf(fidTeX,'\\end{figure}\n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
  	end
	if options_.nograph, close(h1), end
else
	for plt = 1:nbplt-1
		hplt = figure('Name',figurename);
		if TeX
			TeXNAMES = [];
		end    
		NAMES    = []; 
		for index=1:nstar
			names = [];
			i = (plt-1)*nstar + index;
			[borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
			    posterior_distribution(i,nfile,ffil,ifil,...
				     nblck,n,number_of_simulations_per_file,TeX);
	   		eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
	   		eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']);				     
			if TeX
				TeXNAMES = strvcat(TeXNAMES,texnam);
			end    
			NAMES = strvcat(NAMES,nam);
			subplot(nr,nc,index);
			hh = plot(x2,f2,'-k','linewidth',2);
			set(hh,'color',[0.7 0.7 0.7]);
			hold on;
			plot(x1,f1,'-k','linewidth',2);
			plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
			box on;
			axis([borneinf bornesup 0 1.1*top]);
			title(nam,'Interpreter','none');
			hold off;
			drawnow;
		end  % index=1:nstar
		eval(['print -depsc2 ' M_.fname '_PriorsAndPosteriors' int2str(plt)]);
		eval(['print -dpdf ' M_.fname '_PriorsAndPosteriors' int2str(plt)]);
		saveas(hplt,[M_.fname '_PriorsAndPosteriors' int2str(plt) '.fig']);
		if TeX
			fprintf(fidTeX,'\\begin{figure}[H]\n');
			for jj = 1:nstar
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
			end    
			fprintf(fidTeX,'\\centering\n');
			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',M_.fname,int2str(plt));
			fprintf(fidTeX,'\\caption{Priors and posteriors.}');
			fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(plt));
			fprintf(fidTeX,'\\end{figure}\n');
			fprintf(fidTeX,' \n');
		end    
		if options_.nograph, close(hplt), end
	end % plt = 1:nbplt-1
  	hplt = figure('Name',figurename);
  	if TeX
		TeXNAMES = [];
  	end    
  	NAMES    = []; 
  	for index=1:npar-(nbplt-1)*nstar
		i = (nbplt-1)*nstar +  index;
		[borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
    		posterior_distribution(i,nfile,ffil,ifil,...
			   nblck,n,number_of_simulations_per_file,TeX);
	   	eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
	   	eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']);			   
		if TeX
  			TeXNAMES = strvcat(TeXNAMES,texnam);
		end
		NAMES = strvcat(NAMES,nam);
		if lr
  			subplot(lc,lr,index);
		else
  			subplot(nr,nc,index);
		end    
		hh = plot(x2,f2,'-k','linewidth',2);
		set(hh,'color',[0.7 0.7 0.7]);
		hold on;
		plot(x1,f1,'-k','linewidth',2);
		plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
		box on;
		axis([borneinf bornesup 0 1.1*top]);
		title(nam,'Interpreter','none');
		hold off;
		drawnow;
  	end  % index=1:npar-(nbplt-1)*nstar
  	eval(['print -depsc2 ' M_.fname '_PriorsAndPosteriors' int2str(nbplt)]);
  	eval(['print -dpdf ' M_.fname '_PriorsAndPosteriors' int2str(nbplt)]);
  	saveas(hplt,[M_.fname '_PriorsAndPosteriors' int2str(nbplt) '.fig']);
  	if TeX
		fprintf(fidTeX,'\\begin{figure}[H]\n');
		for jj = 1:npar-(nbplt-1)*nstar
			fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
		end    
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',M_.fname,int2str(nbplt));
		fprintf(fidTeX,'\\caption{Priors and posteriors.}');
		fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(nbplt));
		fprintf(fidTeX,'\\end{figure}\n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'%% End of TeX file.\n');
		fclose(fidTeX);
  	end
  	if options_.nograph, close(hplt), end
end
%
%%
%%%
%%%%
%%%%% Je (re)fais mes comptes... I should be able to skip this part (already done) 
%%%%                                                            
%%%                                                             
%%                                                              
%                                                               
FLN = zeros(nfile-ffil+1,3);% Describes the number of lines in each file
if nblck == 1
	instr1 = [M_.fname '_mh'];
	instr2 = '';
else	% I only consider draws from the first chain. This is correct
		% if and only if the metropolis-hastings did converge.
	instr1 = [M_.fname '_mh'];
	instr2 = '_blck1';
end	
eval(['load ' instr1 int2str(ffil) instr2]);
clear post2 x2;
FLN(1,1) = ffil;                            % File number   
FLN(1,2) = size(logpo2(ifil:end,:),1);      % Number of simulations in this file (density) 
FLN(1,3) = FLN(1,2);                        % Cumulative Distribution Function
if nfile-ffil+1>1
	linee = 1;
  	for n = ffil+1:nfile
		linee = linee+1;
		instr = [instr1 int2str(n) instr2];
		eval(['load ' instr]);
		clear post2 x2;
		FLN(linee,1) = n;
		FLN(linee,2) = size(logpo2,1);
		FLN(linee,3) = FLN(linee-1,3) + FLN(linee,2);  
  	end
  	clear logpo2
  	nruns = FLN(linee,3);
else
	nruns = FLN(1,3);    
end
FLN(:,3) = FLN(:,3)/nruns;% I'm scaling the CDF
nvar     = length(oo_.steady_state);
B        = round(0.25*nruns);
deciles = [round(0.1*B) ...
       round(0.2*B)...
       round(0.3*B)...
       round(0.4*B)...
       round(0.5*B)...
       round(0.6*B)...
       round(0.7*B)...
       round(0.8*B)...
       round(0.9*B)];
%                                                               
%%                                                              
%%%                                                             
%%%%                                                            
%%%%% SDGE-based forecasts, smooth and filtered variables, IRFs and theoretical moments 
%%%%                                                            
%%%                                                             
%%                                                              
%                                                               
if options_.forecast | options_.smoother | options_.filtered_vars
    % [1] I delete some old files...    
	disp(' ')
    disp(' ')
    if options_.forecast
    	files = eval(['dir(''' M_.fname '_forecast*.mat'');']);
    	if length(files)
			delete([M_.fname '_forecast*.mat']);
			disp(['MH: Old ' M_.fname '_forecast files deleted! '])
    	end
    end
    if options_.smoother 		
    	files = eval(['dir(''' M_.fname '_smooth*.mat'');']);
    	if length(files)
			delete([M_.fname '_smooth*.mat']);
			disp(['MH: Old ' M_.fname '_smooth files deleted! '])
    	end
    	files = eval(['dir(''' M_.fname '_innovation*.mat'');']);
    	if length(files)
			delete([M_.fname '_innovation*.mat']);
			disp(['MH: Old ' M_.fname '_innovation files deleted! '])
    	end
    	files = eval(['dir(''' M_.fname '_error*.mat'');']);
    	if length(files)
			delete([M_.fname '_error*.mat']);
			disp(['MH: Old ' M_.fname '_error files deleted! '])
    	end
    end
    if options_.filtered_vars
    	files = eval(['dir(''' M_.fname '_filter*.mat'');']);     
    	if length(files)                                        
			delete([M_.fname '_filter*.mat']);                    
			disp(['MH: Old ' M_.fname '_filter files deleted! ']) 
    	end                                                         
    end
    disp(' ')
    disp(' ')
    % [2] Initialization...    
    oo_.exo_simul      = zeros(horizon+M_.maximum_lag+M_.maximum_lead,M_.exo_nbr);
    yyyy     = zeros(nvar,M_.maximum_lag);
    IdObs    = zeros(nvobs,1);
    if options_.forecast 
    	if B <= MAX_nforc
			stock_forcst = zeros(options_.forecast+M_.maximum_lag,nvar,B);
    	else
			stock_forcst = zeros(options_.forecast+M_.maximum_lag,nvar,MAX_nforc);
    	end
    end	
    if options_.smoother
    	if B <= MAX_nsmoo
			stock_smooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,B);
    	else
			stock_smooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,MAX_nsmoo);
    	end
    	if B <= MAX_ninno	
			stock_innov  = zeros(M_.exo_nbr,gend,B);
    	else
			stock_innov  = zeros(M_.exo_nbr,gend,MAX_ninno);
    	end
    	if nvn & B <= MAX_nerro
			stock_error = zeros(gend,nvobs,B);
    	else nvn & B > MAX_nerro
			stock_error = zeros(gend,nvobs,MAX_nerro);
    	end
    end
    if options_.filtered_vars
    	if B <= MAX_nfilt
			stock_filter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,B);
    	else
			stock_filter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,MAX_nfilt);
    	end
    end
    for j=1:nvobs
		for i=1:nvar
	  		iobs = strmatch(options_.varobs(j,:),M_.endo_names,'exact');
		end
		IdObs(j,1) = iobs; 
    end    
    h = waitbar(0,'SDGE model based forecasts...');
    % [3] 	CoRe    
	% [3.1]	First we consider the case with measurement error
    if nvn
		% [3.1.1] More than one _mh file 
		if nfile-ffil+1>1			
  			if options_.forecast
  				sfil_forc = 1;
  				irun_forc = 0;  			
  			end
  			if options_.smoother
  				sfil_smoo = 1;
  				sfil_inno = 1;
  				sfil_erro = 1;
  				irun_smoo = 0;
  				irun_inno = 0;
  				irun_erro = 0;
  			end
  			if options_.filtered_vars
  				sfil_filt = 1;
  				irun_filt = 0;  			
  			end
  			% [3.1.1.1] Loop in the metropolis
  			for b = 1:B;
	    		if options_.forecast
	    			irun_forc = irun_forc+1;
	    		end
	    		if options_.smoother
	    			irun_smoo = irun_smoo+1;
	    			irun_inno = irun_inno+1;
	    			irun_erro = irun_erro+1;
	    		end
  				if options_.filtered_vars
  					irun_filt = irun_filt+1;  			
  				end    			
  				% FIRST, I choose an _mh file (where the posterior distribution is stored)
    			choose_an_mh_file = rand;
    			mh_file_number = FLN(find(choose_an_mh_file>=FLN(:,3)),1);
    			if isempty(mh_file_number)
      				mh_file_number = ffil;
    			else    
      				mh_file_number = mh_file_number(1);
    			end    
    			eval(['load ' instr1 int2str(mh_file_number) instr2]);
    			clear post2 logpo2;
    			% SECOND, I choose a vector of structural parameters (a line in the _mh file) 
    			deep  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
    			% THIRD, I estimate the smooth and filtered variables. I need the smoothed variables
    			% to estimate the state of the model at the end of the sample. 
    			[atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
    			% FOURTH, smoothed and filtered variables are saved if needed
    			if options_.smoother
    				if irun_erro < MAX_nerro
      					stock_error(:,:,irun_erro) = transpose(obs_err);
    				else
      					stock_error(:,:,irun_erro) = transpose(obs_err);
      					instr = [M_.fname '_error' int2str(sfil_erro) ' stock_error;'];
      					eval(['save ' instr]);
      					sfil_erro = sfil_erro + 1;
      					irun_erro = 0;
      					stock_error  = zeros(gend,nvobs,MAX_nerro);
    				end
    				if irun_smoo < MAX_nsmoo
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
    				else
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
      					instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
      					eval(['save ' instr]);
      					sfil_smoo = sfil_smoo + 1;
      					irun_smoo = 0;
      					stock_smooth = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,MAX_nsmoo);
    				end	
    				if irun_inno < MAX_ninno
      					stock_innov(:,:,irun_inno) = innov;
    				else
      					stock_innov(:,:,irun_inno) = innov;
      					instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
      					eval(['save ' instr]);
      					sfil_inno = sfil_inno + 1;
      					irun_inno = 0;
      					stock_innov  = zeros(M_.exo_nbr,gend,MAX_ninno);
    				end	
    			end
    			if options_.filtered_vars
    				if irun_filt < MAX_nfilt
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
    				else
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
      					instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
      					eval(['save ' instr]);
      					sfil_filt = sfil_filt + 1;
      					irun_filt = 0;
      					stock_filter = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,MAX_nfilt);
    				end	
    			end    			
    			if options_.forecast
    				% FIFTH, I update variable dr 
    				dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
    				% SIXTH, I do and save the forecasts (for all the endogenous variables)
    				for j = 1:nvar	% The state of the economy at the end of the sample 
									% depends on the structural parameters.
						if any(j==IdObs)
			  				idx = find(j==IdObs);
			  				yyyy(dr.order_var(j),1:M_.maximum_lag) = data(idx,size(data,2)-M_.maximum_lag+1:end);
						else
			  				yyyy(dr.order_var(j),1:M_.maximum_lag) = atT(j,size(atT,2)-M_.maximum_lag+1:size(atT,2));
						end    
    				end
    				if irun_forc < MAX_nforc
      					stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
    				else
      					stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
      					instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
      					eval(['save ' instr]);
      					sfil_forc = sfil_forc + 1;
      					irun_forc = 0;
      					stock_forcst = zeros(horizon+M_.maximum_lag,nvar,MAX_nforc);
    				end
    			end		
    			waitbar(b/B,h);    
  			end % of loop [3.1.1.1]
  			if options_.smoother
  				if irun_smoo
  				  	stock_smooth = stock_smooth(:,:,1:irun_smoo);
  				  	instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_smooth;
  				if irun_inno
  				  	stock_innov = stock_innov(:,:,1:irun_inno);
  				  	instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_innov;
  				if irun_erro
  				  stock_error = stock_error(:,:,1:irun_erro);
  				  instr = [M_.fname '_error' int2str(sfil_erro) ' stock_error;'];
  				  eval(['save ' instr]);
  				end
  				clear stock_error;
  			end
  			if options_.forecast	
  				if irun_forc
  				  	stock_forcst = stock_forcst(:,:,1:irun_forc);  
  				  	instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_forcst;
  			end
  			if options_.filtered_vars	
  				if irun_filt
  				  stock_filter = stock_filter(:,:,1:irun_filt);
  				  instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
  				  eval(['save ' instr]);
  				end
  				clear stock_filter;
  			end
		else % [3.1.2] Just one _mh file
  			if options_.forecast
  				sfil_forc = 1;
  				irun_forc = 0;  			
  			end
  			if options_.smoother
  				sfil_smoo = 1;
  				sfil_inno = 1;
  				sfil_erro = 1;
  				irun_smoo = 0;
  				irun_inno = 0;
  				irun_erro = 0;
  			end
  			if options_.filtered_vars
  				sfil_filt = 1;
  				irun_filt = 0;  			
  			end
	  		eval(['load ' instr1 int2str(ffil) instr2]);
	  		NumberOfSimulations = length(logpo2);
	  		clear post2 logpo2;
	  		for b = 1:B;
	    		if options_.forecast
	    			irun_forc = irun_forc+1;
	    		end
	    		if options_.smoother
	    			irun_smoo = irun_smoo+1;
	    			irun_inno = irun_inno+1;
	    			irun_erro = irun_erro+1;
	    		end
  				if options_.filtered_vars
  					irun_filt = irun_filt+1;  			
  				end
	    		deep  = x2(floor(rand*NumberOfSimulations)+1,:); 
	    		[atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
	    		if options_.smoother
	    			if irun_erro < MAX_nerro
	    				stock_error(:,:,irun_erro) = transpose(obs_err);
	    			else
	    			  	stock_error(:,:,irun_erro) = transpose(obs_err);
	    			  	instr = [M_.fname '_error' int2str(sfil_erro) ' stock_error;'];
	    			  	eval(['save ' instr]);
	    			  	sfil_erro = sfil_erro + 1;
	    			  	irun_erro = 0;
	    			  	stock_error  = zeros(gend,nvobs,MAX_nerro);
	    			end
	    			if irun_smoo < MAX_nsmoo
	    			  	stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
	    			else
	    			  	stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
	    			  	instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    			  	eval(['save ' instr]);
	    			  	sfil_smoo = sfil_smoo + 1;
	    			  	irun_smoo = 0;
	    			  	stock_smooth = ...
					  	zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,MAX_nsmoo);
	    			end	
	    			if irun_inno < MAX_ninno
	    			  	stock_innov(:,:,irun_inno) = innov;
	    			else
	    			  	stock_innov(:,:,irun_inno) = innov;
	    			  	instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    			  	eval(['save ' instr]);
	    			  	sfil_inno = sfil_inno + 1;
	    			  	irun_inno = 0;
	    			  	stock_innov  = zeros(M_.exo_nbr,gend,MAX_ninno);
	    			end
	    		end
	    		if options_.filtered_vars
    				if irun_filt < MAX_nfilt                                             
      					stock_filter(:,:,irun_filt) = filtered_state_vector;             
    				else                                                                 
      					stock_filter(:,:,irun_filt) = filtered_state_vector;             
      					instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];  
      					eval(['save ' instr]);                                           
      					sfil_filt = sfil_filt + 1;                                       
      					irun_filt = 0;                                                   
      					stock_filter = ...                                               
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,MAX_nfilt);     
    				end	                                                                 
	    		end
	    		if options_.forecast
	    			dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
	    			for j = 1:nvar 
						if any(j==IdObs)
							idx = find(j==IdObs);
							yyyy(dr.order_var(j),1:M_.maximum_lag) = data(idx,size(data,2)-M_.maximum_lag+1:end);
						else
							yyyy(dr.order_var(j),1:M_.maximum_lag) = atT(j,size(atT,2)-M_.maximum_lag+1:size(atT,2));
						end
					end	    
	    			if irun_forc < MAX_nforc
	    		  		stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	    			else
	    		  		stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	    		  		instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
	    		  		eval(['save ' instr]);
	    		  		sfil_forc = sfil_forc + 1;
	    		  		irun_forc = 0;
	    		  		stock_forcst = zeros(horizon+M_.maximum_lag,nvar,MAX_nforc);
	    			end
	    		end	
	    		waitbar(b/B,h);    
	  		end % of the loop over the metropolis simulations
	  		if options_.smoother
	    		if irun_smoo
	    			stock_smooth = stock_smooth(:,:,1:irun_smoo);
	   			 	instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    			eval(['save ' instr]);
	  			end
	  			clear stock_smooth;
	  			if irun_inno
	    			stock_innov = stock_innov(:,:,1:irun_inno);
	    			instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    			eval(['save ' instr]);
	  			end
	  			clear stock_innov;
	  			if irun_erro
	    			stock_error = stock_error(:,:,1:irun_erro);
	    			instr = [M_.fname '_error' int2str(sfil_erro) ' stock_error;'];
	    			eval(['save ' instr]);
	  			end
	  			clear stock_error;
	  		end
	  		if options_.forecast	
	  			if irun_forc
	    			stock_forcst = stock_forcst(:,:,1:irun_forc);  
	    			instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
	    			eval(['save ' instr]);
	  			end
	  			clear stock_forcst;
	  		end
	  		if options_.filtered_vars
  				if irun_filt
  					stock_filter = stock_filter(:,:,1:irun_filt);
  				  	instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_filter;	  		
	  		end
		end
    else % [3.2]	Second we consider the case without measurement error
		if nfile-ffil+1>1
  			if options_.forecast
  				sfil_forc = 1;
  				irun_forc = 0;  			
  			end
  			if options_.smoother
  				sfil_smoo = 1;
  				sfil_inno = 1;
  				sfil_erro = 1;
  				irun_smoo = 0;
  				irun_inno = 0;
  				irun_erro = 0;
  			end
  			if options_.filtered_vars
  				sfil_filt = 1;
  				irun_filt = 0;  			
  			end
			for b = 1:B;
	    		if options_.forecast
	    			irun_forc = irun_forc+1;
	    		end
	    		if options_.smoother
	    			irun_smoo = irun_smoo+1;
	    			irun_inno = irun_inno+1;
	    			irun_erro = irun_erro+1;
	    		end
  				if options_.filtered_vars
  					irun_filt = irun_filt+1;  			
  				end	    
	    		choose_an_mh_file = rand;
	    		mh_file_number = FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	    		if isempty(mh_file_number)
	      			mh_file_number = ffil;
	    		else    
	      			mh_file_number = mh_file_number(1);
	    		end    
	    		eval(['load ' instr1 int2str(mh_file_number) instr2]);
	    		clear post2 logpo2;
	    		deep  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	    		[atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
    			if options_.smoother
    				%if irun_erro < MAX_nerro
      				%	stock_error(:,:,irun_erro) = obs_err;
    				%else
      				%	stock_error(:,:,irun_erro) = obs_err;
      				%	instr = [M_.fname '_error' int2str(sfil_erro) ' stock_error;'];
      				%	eval(['save ' instr]);
      				%	sfil_erro = sfil_erro + 1;
      				%	irun_erro = 0;
      				%	stock_error  = zeros(gend,nvobs,MAX_nerro);
    				%end
    				if irun_smoo < MAX_nsmoo
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
    				else
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
      					instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
      					eval(['save ' instr]);
      					sfil_smoo = sfil_smoo + 1;
      					irun_smoo = 0;
      					stock_smooth = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,MAX_nsmoo);
    				end	
    				if irun_inno < MAX_ninno
      					stock_innov(:,:,irun_inno) = innov;
    				else
      					stock_innov(:,:,irun_inno) = innov;
      					instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
      					eval(['save ' instr]);
      					sfil_inno = sfil_inno + 1;
      					irun_inno = 0;
      					stock_innov  = zeros(M_.exo_nbr,gend,MAX_ninno);
    				end	
    			end
    			if options_.filtered_vars
    				if irun_filt < MAX_nfilt
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
    				else
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
      					instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
      					eval(['save ' instr]);
      					sfil_filt = sfil_filt + 1;
      					irun_filt = 0;
      					stock_filter = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,MAX_nfilt);
    				end
    			end    			
    			if options_.forecast	    
	    			dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
	    			for j = 1:nvar 
						yyyy(dr.order_var(j),1:M_.maximum_lag) = atT(j,size(atT,2)-M_.maximum_lag+1:size(atT,2));
	    			end
	    			if irun_forc < MAX_nforc
	      				stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	    			else
	      				stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	      				instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
	      				eval(['save ' instr]);
	      				sfil_forc = sfil_forc + 1;
	      				irun_forc = 0;
	      				stock_forcst = zeros(horizon+M_.maximum_lag,nvar,MAX_nforc);
	    			end
	    		end
	    		waitbar(b/B,h);
	  		end
  			if options_.smoother
  				if irun_smoo
  				  	stock_smooth = stock_smooth(:,:,1:irun_smoo);
  				  	instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_smooth;
  				if irun_inno
  				  	stock_innov = stock_innov(:,:,1:irun_inno);
  				  	instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_innov;
  			end
  			if options_.forecast	
  				if irun_forc
  				  	stock_forcst = stock_forcst(:,:,1:irun_forc);  
  				  	instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_forcst;
  			end
  			if options_.filtered_vars
				if irun_filt
  				  	stock_filter = stock_filter(:,:,1:irun_filt);
  				  	instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_filter;
  			end
		else % just one _mh file
  			if options_.forecast
  				sfil_forc = 1;
  				irun_forc = 0;  			
  			end
  			if options_.smoother
  				sfil_smoo = 1;
  				sfil_inno = 1;
  				%sfil_erro = 1;
  				irun_smoo = 0;
  				irun_inno = 0;
  				%irun_erro = 0;
  			end
  			if options_.filtered_vars
  				sfil_filt = 1;
  				irun_filt = 0;  			
  			end	  		
	  		eval(['load ' instr1 int2str(ffil) instr2]);
	  		NumberOfSimulations = length(logpo2);
	  		clear post2 logpo2;
	  		for b = 1:B;
	    		if options_.forecast
	    			irun_forc = irun_forc+1;
	    		end
	    		if options_.smoother
	    			irun_smoo = irun_smoo+1;
	    			irun_inno = irun_inno+1;
	    			%irun_erro = irun_erro+1;
	    		end
  				if options_.filtered_vars
  					irun_filt = irun_filt+1;  			
  				end	    
	    		deep  = x2(floor(rand*NumberOfSimulations)+1,:); 
	    		[atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
    			if options_.smoother
    				if irun_smoo < MAX_nsmoo
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
    				else
      					stock_smooth(:,:,irun_smoo) = atT(:,1:gend);
      					instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
      					eval(['save ' instr]);
      					sfil_smoo = sfil_smoo + 1;
      					irun_smoo = 0;
      					stock_smooth = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,MAX_nsmoo);
    				end	
    				if irun_inno < MAX_ninno
      					stock_innov(:,:,irun_inno) = innov;
    				else
      					stock_innov(:,:,irun_inno) = innov;
      					instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
      					eval(['save ' instr]);
      					sfil_inno = sfil_inno + 1;
      					irun_inno = 0;
      					stock_innov  = zeros(M_.exo_nbr,gend,MAX_ninno);
    				end	
    			end
    			if options_.filtered_vars
    				if irun_filt < MAX_nfilt
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
    				else
      					stock_filter(:,:,irun_filt) = filtered_state_vector;
      					instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
      					eval(['save ' instr]);
      					sfil_filt = sfil_filt + 1;
      					irun_filt = 0;
      					stock_filter = ...
	  					zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend+1,MAX_nfilt);
    				end	
    			end    			
    			if options_.forecast	    
	    			dr = resol(oo_.steady_state,options_.oo_.dralgo,options_.linear,options_.order);
	    			for j = 1:nvar 
						yyyy(dr.order_var(j),1:M_.maximum_lag) = atT(j,size(atT,2)-M_.maximum_lag+1:size(atT,2));
	    			end
	    			if irun_forc < MAX_nforc
	      				stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	    			else
	      				stock_forcst(:,:,irun_forc) = transpose(simult_(yyyy,dr,oo_.exo_simul,options_.order));
	      				instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
	      				eval(['save ' instr]);
	      				sfil_forc = sfil_forc + 1;
	      				irun_forc = 0;
	      				stock_forcst = zeros(horizon+M_.maximum_lag,nvar,MAX_nforc);
	    			end
	    		end	
	    		waitbar(b/B,h);    
	  		end
  			if options_.smoother
  				if irun_smoo
  				  	stock_smooth = stock_smooth(:,:,1:irun_smoo);
  				  	instr = [M_.fname '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_smooth;
  				if irun_inno
  				  	stock_innov = stock_innov(:,:,1:irun_inno);
  				  	instr = [M_.fname '_innovation' int2str(sfil_inno) ' stock_innov;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_innov;
  			end
  			if options_.forecast	
  				if irun_forc
  				  	stock_forcst = stock_forcst(:,:,1:irun_forc);  
  				  	instr = [M_.fname '_forecast' int2str(sfil_forc) ' stock_forcst;'];
  				  	eval(['save ' instr]);
  				end
  				clear stock_forcst;
  			end
  			if options_.filtered_vars	
  				if irun_filt
  				  stock_filter = stock_filter(:,:,1:irun_filt);
  				  instr = [M_.fname '_filter' int2str(sfil_filt) ' stock_filter;'];
  				  eval(['save ' instr]);
  				end
  				clear stock_filter;
  			end
		end
    end
    close(h);
end
%%
%% Only a subset of variables may be treated    
%%
varlist = options_.varlist;
if isempty(varlist)
	varlist = M_.endo_names;
	nvar	= size(M_.endo_names,1);
	SelecVariables = transpose(1:nvar);
else
	nvar = size(varlist,1);
	SelecVariables = [];
	for i=1:nvar
		if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
			SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
		end
	end
	IdObs    = zeros(nvobs,1);
    for j=1:nvobs
		for i=1:nvar
	  		iobs = strmatch(options_.varobs(j,:),varlist,'exact');
		end
		if ~isempty(iobs)
			IdObs(j,1) = iobs;
		end
    end
end
if TeX
	varlist_TeX = [];
	for i=1:nvar
		varlist_TeX = strvcat(varlist_TeX,M_.endo_names_tex(SelecVariables(i),:));
	end
end
%%                                    %%
%% Now I treat the forecasts (plots)  %%   
%%                                    %%
if options_.forecast
	tmp = zeros(B,1);
	fprintf('MH: Out of sample forecasts...\n');
  	MeanForecast = zeros(options_.forecast,nvar);
  	MedianForecast = zeros(options_.forecast,nvar);
  	StdForecast = zeros(options_.forecast,nvar);
  	HPD   = zeros(options_.forecast,nvar,2);
  	for step = 1:options_.forecast % ... Suffering is one very long moment.
		truestep = step+M_.maximum_lag;
		for i = 1:nvar;
  			StartLine = 0;
  			for file = 1:sfil_forc;
    			instr = [M_.fname '_forecast' int2str(file)];
    			eval(['load ' instr]);
    			MeanForecast(step,i) = MeanForecast(step,i)+sum(stock_forcst(truestep,SelecVariables(i),:),3);
    			DeProfundis = size(stock_forcst,3); 
    			tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_forcst(truestep,SelecVariables(i),:)); 
    			StartLine = StartLine+DeProfundis;
  			end
  			tmp = sort(tmp);
  			MedianForecast(step,i) = tmp(round(B*0.5));
  			StdForecast(step,i) = std(tmp);
  			t = floor(options_.mh_conf_sig*B);
  			a = 1; 
  			b = t;
  			tmp2 = [1;t;tmp(t)-tmp(1)];
  			while b <= B
    			tmp1 = [a;b;tmp(b)-tmp(a)];
    			a = a + 1;
    			b = b + 1;
    			if tmp1(3,1) < tmp2(3,1)
      				tmp2 = tmp1;     
    			end    
  			end
  			HPD(step,i,1) = tmp(tmp2(1,1));
  			HPD(step,i,2) = tmp(tmp2(2,1));
		end
		disp(['    Period = ' int2str(step)]);
  	end
  	MeanForecast = MeanForecast/B;
  	[nbplt,nr,nc,lr,lc,nstar] = pltorg(nvar);
  	if TeX
		fidTeX = fopen([M_.fname '_BayesianForecasts.TeX'],'w');
		fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
		fprintf(fidTeX,' \n');
		NAMES = [];
		TEXNAMES = [];
  	end    
  	if nbplt == 1
		hfig = figure('Name','Out of sample forecasts');
		for i = 1:nvar
  			subplot(nr,nc,i)
  			if any(i==IdObs)
    			idx = find(i==IdObs);
    			plot(1:options_.forecast+10,[data(idx,size(data,2)-10+1:end)';MeanForecast(:,i)],'-b','linewidth',2)
    			hold on            
    			plot(11:10+options_.forecast,HPD(:,i,1),'--k','linewidth',1.5)
    			plot(11:10+options_.forecast,HPD(:,i,2),'--k','linewidth',1.5)
    			plot([11 11],ylim,'-g')
    			set(gca,'XTick',[11 20 30 40 50 60 70 80 90 100]);
    			set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
    			xlim([1 options_.forecast+10]);
    			box on
    			hold off
    			title(deblank(varlist(i,:)),'Interpreter','none')
    			eval(['oo_.Forecast.Mean.' deblank(varlist(i,:)) ' = MeanForecast(:,i)'';']);
    			eval(['oo_.Forecast.Median.' deblank(varlist(i,:)) ' = MedianForecast(:,i)'';']);
    			eval(['oo_.Forecast.Std.' deblank(varlist(i,:)) ' = StdForecast(:,i)'';']);
    			eval(['oo_.Forecast.HPDinf.' deblank(varlist(i,:)) ' = squeeze(HPD(:,i,1))'';']);
    			eval(['oo_.Forecast.HPDsup.' deblank(varlist(i,:)) ' = squeeze(HPD(:,i,2))'';']);
    			if TeX
      				NAMES = strvcat(NAMES,deblank(varlist(i,:)));
      				TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(i,:)) ' $']);
    			end
  			else    
    			plot(1:options_.forecast,HPD(:,i,1),'--k','linewidth',1.5)
    			hold on
    			plot(1:options_.forecast,HPD(:,i,2),'--k','linewidth',1.5)
    			plot(1:options_.forecast,MeanForecast(:,i),'-b','linewidth',2)
    			set(gca,'XTick',[1 10 20 30 40 50 60 70 80 90]);
    			set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});    			
    			xlim([1 options_.forecast]);
    			box on
    			hold off
    			title(deblank(varlist(i,:)),'Interpreter','none')
    			eval(['oo_.Forecast.Mean.' deblank(varlist(i,:)) ' = MeanForecast(:,i)'';']);
    			eval(['oo_.Forecast.Median.' deblank(varlist(i,:)) ' = MedianForecast(:,i)'';']);
    			eval(['oo_.Forecast.Std.' deblank(varlist(i,:)) ' = StdForecast(:,i)'';']);
    			eval(['oo_.Forecast.HPDinf.' deblank(varlist(i,:)) ' = squeeze(HPD(:,i,1))'';']);
    			eval(['oo_.Forecast.HPDsup.' deblank(varlist(i,:)) ' = squeeze(HPD(:,i,2))'';']);
    			if TeX
    				NAMES = strvcat(NAMES,deblank(varlist(i,:)));
    				TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(i,:)) ' $']);
    			end
  			end   
		end
		eval(['print -depsc2 ' M_.fname '_Forecasts' int2str(1)]);
		eval(['print -dpdf ' M_.fname '_Forecasts' int2str(1)]);
		saveas(hfig,[M_.fname '_Forecasts' int2str(1) '.fig']);
		if options_.nograph, close(hfig), end
		if TeX
  			fprintf(fidTeX,'\\begin{figure}[H]\n');
  			for jj = 1:nvar
    			fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
  			end    
  			fprintf(fidTeX,'\\centering \n');
  			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Forecasts%s}\n',M_.fname,int2str(1));
  			fprintf(fidTeX,'\\caption{DSGE posterior mean forecats with HPD intervals.}');
  			fprintf(fidTeX,'\\label{Fig:Forecasts:%s}\n',int2str(1));
  			fprintf(fidTeX,'\\end{figure}\n');
  			fprintf(fidTeX,' \n');
  			fprintf(fidTeX,' \n');
  			fprintf(fidTeX,'% End Of TeX File. \n');
  			fclose(fidTeX); 
		end
	else
		for plt = 1:nbplt-1
  			if TeX
    			NAMES = [];
    			TEXNAMES = [];
  			end
  			hfig = figure('Name','Out of sample forecasts');
  			for i = 1:nstar
    			k = (plt-1)*nstar+i;
    			subplot(nr,nc,i)
    			if any(k==IdObs)
      				idx = find(k==IdObs);
      				plot(1:options_.forecast+10,[data(idx,size(data,2)-10+1:end)';MeanForecast(:,k)],'-b','linewidth',2)
      				hold on            
      				plot(11:10+options_.forecast,HPD(:,k,1),'--k','linewidth',1.5)
      				plot(11:10+options_.forecast,HPD(:,k,2),'--k','linewidth',1.5)
      				plot([11 11],ylim,'-g')
    				set(gca,'XTick',[11 20 30 40 50 60 70 80 90 100]);
    				set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
    				xlim([1 options_.forecast+10]);
      				box on
      				title(deblank(varlist(k,:)),'Interpreter','none')
      				hold off
    				eval(['oo_.Forecast.Mean.' deblank(varlist(k,:)) ' = MeanForecast(:,k)'';']);
    				eval(['oo_.Forecast.Median.' deblank(varlist(k,:)) ' = MedianForecast(:,k)'';']);
    				eval(['oo_.Forecast.Std.' deblank(varlist(k,:)) ' = StdForecast(:,k)'';']);
    				eval(['oo_.Forecast.HPDinf.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,1))'';']);
    				eval(['oo_.Forecast.HPDsup.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,2))'';']);
      				if TeX
						NAMES = strvcat(NAMES,deblank(varlist(k,:)));
						TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(k,:)) ' $']);
      				end
    			else    
      				plot(1:options_.forecast,HPD(:,k,1),'--k','linewidth',1.5)
      				hold on
      				plot(1:options_.forecast,HPD(:,k,2),'--k','linewidth',1.5)
      				plot(1:options_.forecast,MeanForecast(:,k),'-b','linewidth',2)
    				set(gca,'XTick',[1 10 20 30 40 50 60 70 80 90]);
    				set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});    			
    				xlim([1 options_.forecast]);
      				box on
      				title(deblank(varlist(k,:)),'Interpreter','none')
      				hold off
    				eval(['oo_.Forecast.Mean.' deblank(varlist(k,:)) ' = MeanForecast(:,k)'';']);
    				eval(['oo_.Forecast.Median.' deblank(varlist(k,:)) ' = MedianForecast(:,k)'';']);
    				eval(['oo_.Forecast.Std.' deblank(varlist(k,:)) ' = StdForecast(:,k)'';']);
    				eval(['oo_.Forecast.HPDinf.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,1))'';']);
    				eval(['oo_.Forecast.HPDsup.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,2))'';']);
      				if TeX
						NAMES = strvcat(NAMES,deblank(varlist(k,:)));
						TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(k,:)) ' $']);
      				end
    			end   
  			end
  			eval(['print -depsc2 ' M_.fname '_Forecasts' int2str(plt)]);
  			eval(['print -dpdf ' M_.fname '_Forecasts' int2str(plt)]);
  			saveas(hfig,[M_.fname '_Forecasts' int2str(plt) '.fig']);
  			if options_.nograph, close(hfig), end
  			if TeX
    			fprintf(fidTeX,'\\begin{figure}[H]\n');
    			for jj = 1:nstar
      				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
    			end    
    			fprintf(fidTeX,'\\centering \n');
    			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Forecasts%s}\n',M_.fname,int2str(plt));
    			fprintf(fidTeX,'\\caption{DSGE posterior mean forecats with HPD intervals.}');
    			fprintf(fidTeX,'\\label{Fig:Forecasts:%s}\n',int2str(plt));
    			fprintf(fidTeX,'\\end{figure}\n');
    			fprintf(fidTeX,' \n');
  			end
		end
		hfig = figure('Name','Out of sample forecasts');
		if TeX
  			NAMES = [];
  			TEXNAMES = [];
		end
		for i=1:nfor-(nbplt-1)*nstar
  			k = (nbplt-1)*nstar+i;
  			subplot(lr,lc,i);
  			if any(k==IdObs)
    			idx = find(k==IdObs);
    			plot(1:options_.forecast+10,[data(idx,size(data,2)-10+1:end)';MeanForecast(:,k)],'-b','linewidth',2)
    			hold on            
    			plot(11:10+options_.forecast,HPD(:,k,1),'--k','linewidth',1.5)
    			plot(11:10+options_.forecast,HPD(:,k,2),'--k','linewidth',1.5)
    			plot([11 11],'-g')
    			set(gca,'XTick',[11 20 30 40 50 60 70 80 90 100]);
    			set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
    			xlim([1 options_.forecast+10]);
    			box on
    			title(deblank(varlist(k,:)),'Interpreter','none')
    			hold off
    			eval(['oo_.Forecast.Mean.' deblank(varlist(k,:)) ' = MeanForecast(:,k)'';']);
    			eval(['oo_.Forecast.Median.' deblank(varlist(k,:)) ' = MedianForecast(:,k)'';']);
    			eval(['oo_.Forecast.Std.' deblank(varlist(k,:)) ' = StdForecast(:,k)'';']);
    			eval(['oo_.Forecast.HPDinf.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,1))'';']);
    			eval(['oo_.Forecast.HPDsup.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,2))'';']);
    			if TeX
      				NAMES = strvcat(NAMES,deblank(varlist(k,:)));
      				TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(k,:)) ' $']);
      			end	
    		else
				plot(1:options_.forecast,HPD(:,k,1),'--k','linewidth',1.5)                               
    			hold on                                                                                  
    			plot(1:options_.forecast,HPD(:,k,2),'--k','linewidth',1.5)                               
    			plot(1:options_.forecast,MeanForecast(:,k),'-b','linewidth',2)                           
				set(gca,'XTick',[1 10 20 30 40 50 60 70 80 90]);                                         
				set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});    			 
				xlim([1 options_.forecast]);                                                             
    			box on                                                                                   
    			title(deblank(varlist(k,:)),'Interpreter','none')                                        
    			hold off                                                                                 
    			eval(['oo_.Forecast.Mean.' deblank(varlist(k,:)) ' = MeanForecast(:,k)'';']);
    			eval(['oo_.Forecast.Median.' deblank(varlist(k,:)) ' = MedianForecast(:,k)'';']);
    			eval(['oo_.Forecast.Std.' deblank(varlist(k,:)) ' = StdForecast(:,k)'';']);
    			eval(['oo_.Forecast.HPDinf.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,1))'';']);
    			eval(['oo_.Forecast.HPDsup.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,2))'';']);    			
    			if TeX                                                                               
      				NAMES = strvcat(NAMES,deblank(varlist(k,:)));                                    
      				TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(k,:)) ' $']);              
    			end                                                                                  
    		end
    	end	 
		eval(['print -depsc2 ' M_.fname '_Forecasts' int2str(nbplt)]);
		eval(['print -dpdf ' M_.fname '_Forecasts' int2str(nbplt)]);
		saveas(hfig,[M_.fname '_Forecasts' int2str(nbplt) '.fig']);
		if options_.nograph, close(hfig), end
		if TeX
  			fprintf(fidTeX,'\\begin{figure}[H]\n');
  			for jj = 1:nstar
    			fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
  			end    
  			fprintf(fidTeX,'\\centering \n');
  			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Forecasts%s}\n',M_.fname,int2str(nplt));
  			fprintf(fidTeX,'\\caption{DSGE posterior mean forecats with HPD intervals.}');
 	 		fprintf(fidTeX,'\\label{Fig:Forecasts:%s}\n',int2str(nplt));
  			fprintf(fidTeX,'\\end{figure}\n');
  			fprintf(fidTeX,' \n');
  			fprintf(fidTeX,' \n');
  			fprintf(fidTeX,'% End Of TeX File. \n');
  			fclose(fidTeX);
		end
	end
	fprintf('MH: Out of sample forecasts, done!\n')
	disp(' ')
end
%%
%% Smooth variables and Filtered variables (all endogenous variables are considered here)        
%%
if options_.smoother
	fprintf('MH: Smooth variables...\n')
	MeanSmooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	MedianSmooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	StdSmooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	DistribSmooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,9);
	HPDSmooth = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,2);
	for i = 1:size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic;
		for t = 1:gend
			StartLine = 0;
		  	for file = 1:sfil_smoo;
		    	instr = [M_.fname '_smooth' int2str(file)];
		    	eval(['load ' instr]);
		    	MeanSmooth(i,t) = MeanSmooth(i,t)+sum(stock_smooth(i,t,:),3);
		    	DeProfundis = size(stock_smooth,3); 
		    	tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_smooth(i,t,:)); 
		    	StartLine = StartLine+DeProfundis;
		  	end
		  	tmp = sort(tmp);
		  	MedianSmooth(i,t) = tmp(round(B*0.5));
		  	StdSmooth(i,t) = std(tmp);
		  	DistribSmooth(i,t,:) = reshape(tmp(deciles),1,1,9);
		  	tt = floor(options_.mh_conf_sig*B);
		  	a = 1; 
		  	b = tt;
		  	tmp2 = [1;tt;tmp(tt)-tmp(1)];
		  	while b <= B
		    	tmp1 = [a;b;tmp(b)-tmp(a)];
		    	a = a + 1;
		    	b = b + 1;
		    	if tmp1(3,1) < tmp2(3,1)
		      		tmp2 = tmp1;     
		    	end    
		  	end
		  	HPDSmooth(i,t,1) = tmp(tmp2(1,1));
		  	HPDSmooth(i,t,2) = tmp(tmp2(2,1));
		end
		disp(['    Variable: ' deblank(M_.endo_names(oo_.dr.order_var(i),:))]);	
	end
    clear stock_smooth;
    MeanSmooth = MeanSmooth/B;
    for i=1:size(M_.endo_names,1)
		eval(['oo_.PosteriorSmoothedVariables.Mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = MeanSmooth(i,:)'';']);
		eval(['oo_.PosteriorSmoothedVariables.Median.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = MedianSmooth(i,:)'';']);
		eval(['oo_.PosteriorSmoothedVariables.Std.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = StdSmooth(i,:)'';']);
		eval(['oo_.PosteriorSmoothedVariables.Distribution.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(DistribSmooth(i,:,:))'';']);
		eval(['oo_.PosteriorSmoothedVariables.HPDinf.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(HPDSmooth(i,:,1))'';']);
		eval(['oo_.PosteriorSmoothedVariables.HPDsup.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(HPDSmooth(i,:,2))'';']);
	end
    fprintf('MH: Smooth variables, done!\n')
    disp(' ')
    fprintf('MH: Smooth structural shocks...\n')
    MeanInnov = zeros(M_.exo_nbr,gend);
    MedianInnov = zeros(M_.exo_nbr,gend);
    StdInnov = zeros(M_.exo_nbr,gend);
    DistribInnov = zeros(M_.exo_nbr,gend,9);
    HPDInnov = zeros(M_.exo_nbr,gend,2);
    for i = 1:M_.exo_nbr;
		for t = 1:gend
	  		StartLine = 0;
	  		for file = 1:sfil_inno;
	    		instr = [M_.fname '_innovation' int2str(file)];
	    		eval(['load ' instr]);
	    		MeanInnov(i,t) = MeanInnov(i,t)+sum(stock_innov(i,t,:),3);
	    		DeProfundis = size(stock_innov,3); 
	    		tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_innov(i,t,:)); 
	    		StartLine = StartLine+DeProfundis;
	  		end
	  		tmp = sort(tmp);
	  		MedianInnov(i,t) = tmp(round(B*0.5));
	  		StdInnov(i,t) = std(tmp);
	  		DistribInnov(i,t,:) = reshape(tmp(deciles),1,1,9);
	  		tt = floor(options_.mh_conf_sig*B);
	  		a = 1; 
	  		b = tt;
	  		tmp2 = [1;tt;tmp(tt)-tmp(1)];
	  		while b <= B
	    		tmp1 = [a;b;tmp(b)-tmp(a)];
	    		a = a + 1;
	    		b = b + 1;
	    		if tmp1(3,1) < tmp2(3,1)
	      			tmp2 = tmp1;     
	    		end    
	  		end
	  		HPDInnov(i,t,1) = tmp(tmp2(1,1));
	  		HPDInnov(i,t,2) = tmp(tmp2(2,1));
		end
		disp(['    Variable: ' deblank(M_.exo_names(i,:))]);
	end
    clear stock_innov;
    MeanInnov = MeanInnov/B;
    for i=1:M_.exo_nbr
		eval(['oo_.PosteriorSmoothedShocks.Mean.' deblank(M_.exo_names(i,:)) ' = MeanInnov(i,:)'';']);
		eval(['oo_.PosteriorSmoothedShocks.Median.' deblank(M_.exo_names(i,:)) ' = MedianInnov(i,:)'';']);
		eval(['oo_.PosteriorSmoothedShocks.Std.' deblank(M_.exo_names(i,:)) ' = StdInnov(i,:)'';']);
		eval(['oo_.PosteriorSmoothedShocks.Distribution.' deblank(M_.exo_names(i,:)) ' = squeeze(DistribInnov(i,:,:))'';']);
		eval(['oo_.PosteriorSmoothedShocks.HPDinf.' deblank(M_.exo_names(i,:)) ' = squeeze(HPDInnov(i,:,1))'';']);
		eval(['oo_.PosteriorSmoothedShocks.HPDsup.' deblank(M_.exo_names(i,:)) ' = squeeze(HPDInnov(i,:,2))'';']);
    end
    fprintf('MH: Smooth structural shocks, done!\n')
    disp(' ')
    if nvn
		fprintf('MH: Smooth measurement error...\n')
		MeanError = zeros(gend,nvobs);
		MedianError = zeros(gend,nvobs);
		StdError = zeros(gend,nvobs);
		DistribError = zeros(gend,nvobs,9);
		HPDError = zeros(gend,nvobs,2);
		for i = 1:nvobs;
			for t = 1:gend
		    	StartLine = 0;
		    	for file = 1:sfil_erro;
		      		instr = [M_.fname '_error' int2str(file)];
		      		eval(['load ' instr]);
		      		MeanError(t,i) = MeanError(t,i)+sum(stock_error(t,i,:),3);
		      		DeProfundis = size(stock_error,3); 
		      		tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_error(t,i,:)); 
		      		StartLine = StartLine+DeProfundis;
		    	end
		    	tmp = sort(tmp);
		    	MedianError(t,i) = tmp(round(B*0.5));
		    	StdError(t,i) = std(tmp);
		    	DistribError(t,i,:) = reshape(tmp(deciles),1,1,9);
		    	tt = floor(options_.mh_conf_sig*B);
		    	a = 1; 
		    	b = tt;
		    	tmp2 = [1;tt;tmp(tt)-tmp(1)];
		    	while b <= B
		      		tmp1 = [a;b;tmp(b)-tmp(a)];
		      		a = a + 1;
		      		b = b + 1;
		      		if tmp1(3,1) < tmp2(3,1)
						tmp2 = tmp1;     
		      		end    
		    	end
		    	HPDError(t,i,1) = tmp(tmp2(1,1));
		    	HPDError(t,i,2) = tmp(tmp2(2,1));
		  	end
		  	disp(['    Variable: ' deblank(options_.varobs(i))]);
		end
		clear stock_error;
		MeanError = MeanError/B;
		for i=1:nvobs
			eval(['oo_.PosteriorSmoothedMeasurementErrors.Mean.' deblank(options_.varobs(i)) ' = MeanError(:,i)'';']);
		  	eval(['oo_.PosteriorSmoothedMeasurementErrors.Median.' deblank(options_.varobs(i)) ' = MedianError(:,i)'';']);
		  	eval(['oo_.PosteriorSmoothedMeasurementErrors.Std.' deblank(options_.varobs(i)) ' = StdError(:,i)'';']);
		  	eval(['oo_.PosteriorSmoothedMeasurementErrors.Distribution.' deblank(options_.varobs(i)) ' = squeeze(DistribError(:,i,:))'';']);
		  	eval(['oo_.PosteriorSmoothedMeasurementErrors.HPDinf.' deblank(options_.varobs(i)) ' = squeeze(HPDError(:,i,1))'';']);
			eval(['oo_.PosteriorSmoothedMeasurementErrors.HPDsup.' deblank(options_.varobs(i)) ' = squeeze(HPDError(:,i,2))'';']);
		end
		fprintf('MH: Smooth measurement error, done!\n')
		disp(' ')
	end	
    %%
    %% Now I plot the smooth -structural- shocks
    %%
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
    if TeX
		fidTeX = fopen([M_.fname '_SmoothedShocks.TeX'],'w');
		fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
		fprintf(fidTeX,' \n');
	end
    if nbplt == 1
		hh = figure('Name','Smoothed shocks');    
		NAMES = [];
		if TeX; TEXNAMES = []; end;
		for i=1:M_.exo_nbr
	  		set(0,'CurrentFigure',hh)
	  		subplot(nr,nc,i);
	  		plot([1 gend],[0 0],'-r','linewidth',0.5);
	  		hold on
	  		for j = 1:9
	    		plot(1:gend,DistribInnov(i,:,j),'-g','linewidth',0.5)
	  		end
	  		plot(1:gend,MeanInnov(i,:),'-k','linewidth',1)
	  		xlim([1 gend]);
	  		hold off
	  		ih = figure('Visible','off');
	  		set(0,'CurrentFigure',ih)
	  		plot([1 gend],[0 0],'-r','linewidth',0.5);
	  		hold on
	  		for j = 1:9
	    		plot(1:gend,DistribInnov(i,:,j),'-g','linewidth',0.5)
	  		end
	  		plot(1:gend,MeanInnov(i,:),'-k','linewidth',1)
	  		xlim([1 gend]);
	  		hold off
	  		name    = M_.exo_names(i,:);
	  		NAMES   = strvcat(NAMES,name);
	  		if ~isempty(options_.XTick)
	    		set(gca,'XTick',options_.XTick)
	    		set(gca,'XTickLabel',options_.XTickLabel)
	  		end
	  		eval(['print -depsc2 ' M_.fname '_SmoothedShock_' name]);
	  		eval(['print -dpdf ' M_.fname '_SmoothedShock_' name]);
	  		saveas(ih,[M_.fname '_SmoothedShock_' name '.fig']);
	  		if TeX
	    		texname = M_.exo_names_tex(i,1);
	    		TEXNAMES   = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
	  		end
	  		set(0,'CurrentFigure',hh)
	  		title(name,'Interpreter','none')
	  		if ~isempty(options_.XTick)
	    		set(gca,'XTick',options_.XTick)
	    		set(gca,'XTickLabel',options_.XTickLabel)
	  		end
		end
		eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(1)]);
		eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(1)]);
		saveas(hh,[M_.fname '_SmoothedShocks' int2str(1) '.fig']);
		if options_.nograph, close(hh), end
		if TeX
	  		fprintf(fidTeX,'\\begin{figure}[H]\n');
	  		for jj = 1:M_.exo_nbr
	    		fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	  		end    
	  		fprintf(fidTeX,'\\centering \n');
	  		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(1));
	  		fprintf(fidTeX,'\\caption{Smoothed shocks.}');
	  		fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(1));
	  		fprintf(fidTeX,'\\end{figure}\n');
	  		fprintf(fidTeX,' \n');
	  		fprintf(fidTeX,'%% End of TeX file.\n');
	  		fclose(fidTeX);
		end    
	else
		for plt = 1:nbplt-1
	  		hh = figure('Name','Smoothed shocks');
	  		NAMES = [];
	  		if TeX; TEXNAMES = []; end;
	  		for i=1:nstar
	    		k = (plt-1)*nstar+i;
	    		set(0,'CurrentFigure',hh)
	    		subplot(nr,nc,i);
	    		plot([1 gend],[0 0],'-r','linewidth',0.5)
	    		hold on
	    		for j = 1:9
	      			plot(1:gend,DistribInnov(k,:,j),'-g','linewidth',0.5)
	    		end
	    		plot(1:gend,MeanInnov(k,:),'-k','linewidth',1)
	    		xlim([1 gend])
	    		hold off
	    		name = M_.exo_names(k,:);
	    		NAMES = strvcat(NAMES,name);
	    		ih = figure('Visible','off');
	    		set(0,'CurrentFigure',ih)
	    		plot([1 gend],[0 0],'-r','linewidth',0.5);
	    		hold on
	    		for j = 1:9
	    			plot(1:gend,DistribInnov(k,:,j),'-g','linewidth',0.5)
	    		end
	    		plot(1:gend,MeanInnov(k,:),'-k','linewidth',1)
	    		xlim([1 gend]);
	    		hold off
	    		if ~isempty(options_.XTick)
	    			set(gca,'XTick',options_.XTick)
	    			set(gca,'XTickLabel',options_.XTickLabel)
	    		end
	    		eval(['print -depsc2 ' M_.fname '_SmoothedShock_' name]);
	    		eval(['print -dpdf ' M_.fname '_SmoothedShock_' name]);
	    		saveas(ih,[M_.fname '_SmoothedShock_' name '.fig']);
	    		if TeX
	    			texname = M_.exo_names_tex(k,:);
	    			TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
	    		end    
	    		set(0,'CurrentFigure',hh)
	    		title(name,'Interpreter','none')
	    		if ~isempty(options_.XTick)
	    			set(gca,'XTick',options_.XTick)
	    			set(gca,'XTickLabel',options_.XTickLabel)
	    		end
	  		end
	  		eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(plt)]);
	  		eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(plt)]);
	  		saveas(hh,[M_.fname '_SmoothedShocks' int2str(plt) '.fig']);
	  		if options_.nograph, close(hh), end
	  		if TeX
	  			fprintf(fileone,'\\begin{figure}[H]\n');
	  			for jj = 1:nstar
	  			  fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	  			end    
	  			fprintf(fidTeX,'\\centering \n');
	  			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(plt));
	  			fprintf(fidTeX,'\\caption{Smoothed shocks.}');
	  			fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
	  			fprintf(fidTeX,'\\end{figure}\n');
	  			fprintf(fidTeX,' \n');
	  		end    
		end
		hh = figure('Name','Smoothed shocks');
		NAMES = [];
		if TeX; TEXNAMES = []; end;
		for i=1:M_.exo_nbr-(nbplt-1)*nstar
	  		k = (nbplt-1)*nstar+i;
	  		set(0,'CurrentFigure',hh)
	  		if lr ~= 0
	  			subplot(lr,lc,i);
	  		else
	  			subplot(nr,nc,i);
	  		end    
	  		plot([1 gend],[0 0],'-r','linewidth',0.5)
	  		hold on
	  		for j = 1:9
	  			plot(1:gend,DistribInnov(k,:,j),'-g','linewidth',0.5)
	  		end
	  		plot(1:gend,MeanInnov(k,:),'-k','linewidth',1)
	  		xlim([1 gend]);
	  		hold off
	  		name = M_.exo_names(k,:);
	  		NAMES = strvcat(NAMES,name);
	  		ih = figure('Visible','off');
	  		set(0,'CurrentFigure',ih)
	  		plot([1 gend],[0 0],'-r','linewidth',0.5);
	  		hold on
	  		for j = 1:9
	  			plot(1:gend,DistribInnov(k,:,j),'-g','linewidth',0.5)
	  		end
	  		plot(1:gend,MeanInnov(k,:),'-k','linewidth',1)
	  		xlim([1 gend]);
	  		hold off
	  		if ~isempty(options_.XTick)
	  			set(gca,'XTick',options_.XTick)
	  			set(gca,'XTickLabel',options_.XTickLabel)
	  		end
	  		eval(['print -depsc2 ' M_.fname '_SmoothedShock_' name]);
	  		eval(['print -dpdf ' M_.fname '_SmoothedShock_' name]);
	  		saveas(ih,[M_.fname '_SmoothedShock_' name '.fig']);
	  		if TeX
	    		texname  = M_.exo_names_tex(k,:);
	    		TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
	  		end
	  		set(0,'CurrentFigure',hh)
	  		title(name,'Interpreter','none');
	  		if ~isempty(options_.XTick)
	    		set(gca,'XTick',options_.XTick)
	    		set(gca,'XTickLabel',options_.XTickLabel)
	  		end
		end
		eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(nbplt)]);
		eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(nbplt)]);
		saveas(hh,[M_.fname '_SmoothedShocks' int2str(nbplt) '.fig']);
		if options_.nograph, close(hh), end
		if TeX
	  		fprintf(fileone,'\\begin{figure}[H]\n');
	  		for jj = 1:nstar
	    		fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	  		end    
	  		fprintf(fidTeX,'\\centering \n');
	  		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(nbplt));
	  		fprintf(fidTeX,'\\caption{Smoothed shocks.}');
	  		fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(nbplt));
	  		fprintf(fidTeX,'\\end{figure}\n');
	  		fprintf(fidTeX,' \n');
	  		fprintf(fidTeX,' \n');
	  		fprintf(fidTeX,'%% End of TeX file.\n');
	  		fclose(fidTeX);
		end
	end % nbplt == 1 (smooth -structural- shocks)	
	%%
	%%	Smoothed variables (observed and unobserved)
	%%
	%%	Here non zero steady state levels and linear trends are removed... Would be nice
	%%  to add them... Later.
	[nbplt,nr,nc,lr,lc,nstar] = pltorg(size(M_.endo_names,1));
	if TeX
		fidTeX = fopen([M_.fname '_SmoothedVariables.TeX'],'w');
	    fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
	    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
	    fprintf(fidTeX,' \n');
	end    
	if nbplt == 1
	    hh = figure('Name','Smoothed variables');    
	    NAMES = [];
	    if TeX; TEXNAMES = []; end;
	    for i=1:size(M_.endo_names,1)
	    	set(0,'CurrentFigure',hh)
	      	subplot(nr,nc,i);
	      	plot([1 gend],[0 0],'-r','linewidth',0.5);
	      	hold on
	      	for j = 1:9
				plot(1:gend,DistribSmooth(i,:,j),'-g','linewidth',0.5)
	      	end
	      	plot(1:gend,MeanSmooth(i,:),'-k','linewidth',1)
	      	xlim([1 gend]);
	      	hold off
	      	ih = figure('Visible','off');
	      	set(0,'CurrentFigure',ih)
	      	plot([1 gend],[0 0],'-r','linewidth',0.5);
	      	hold on
	      	for j = 1:9
				plot(1:gend,DistribSmooth(i,:,j),'-g','linewidth',0.5)
	      	end
	      	plot(1:gend,MeanSmooth(i,:),'-k','linewidth',1)
	      	xlim([1 gend]);
	      	hold off
	      	name = deblank(M_.endo_names(oo_.dr.order_var(i),:));
	      	NAMES   = strvcat(NAMES,name);
	      	if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
	      	end
	      	eval(['print -depsc2 ' M_.fname '_SmoothedVariable_' name]);
	      	eval(['print -dpdf ' M_.fname '_SmoothedVariable_' name]);
	      	saveas(ih,[M_.fname '_SmoothedVariable_' name '.fig']);
	      	if TeX
				texname = deblank(M_.endo_names_tex(oo_.dr.order_var(i),1));
				TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
	      	end
	      	set(0,'CurrentFigure',hh)
	      	title(name,'Interpreter','none')
	      	if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
	      	end
	    end
	    eval(['print -depsc2 ' M_.fname '_SmoothedVariables' int2str(1)]);
	    eval(['print -dpdf ' M_.fname '_SmoothedVariables' int2str(1)]);
	    saveas(hh,[M_.fname '_SmoothedVariables' int2str(1) '.fig']);
	    if options_.nograph, close(hh), end
	    if TeX
	      	fprintf(fidTeX,'\\begin{figure}[H]\n');
	      	for jj = 1:M_.exo_nbr
				fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
			end    
	      	fprintf(fidTeX,'\\centering \n');
	      	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedVariables%s}\n',M_.fname,int2str(1));
	      	fprintf(fidTeX,'\\caption{Smoothed variables.}');
	      	fprintf(fidTeX,'\\label{Fig:SmoothedVariables:%s}\n',int2str(1));
	      	fprintf(fidTeX,'\\end{figure}\n');
	      	fprintf(fidTeX,' \n');
	      	fprintf(fidTeX,'%% End of TeX file.\n');
	      	fclose(fidTeX);
	    end    
	else
		for plt = 1:nbplt-1
	      	hh = figure('Name','Smoothed variables');
	      	NAMES = [];
	      	TEXNAMES = [];
	      	for i=1:nstar
				k = (plt-1)*nstar+i;
				set(0,'CurrentFigure',hh)
				subplot(nr,nc,i);
				plot([1 gend],[0 0],'-r','linewidth',0.5)
				hold on
				for j = 1:9
					plot(1:gend,DistribSmooth(k,:,j),'-g','linewidth',0.5)
				end
				plot(1:gend,MeanSmooth(k,:),'-k','linewidth',1)
				xlim([1 gend]);
				hold off
				name = deblank(M_.endo_names(oo_.dr.order_var(k),:));
				NAMES = strvcat(NAMES,name);
				ih = figure('Visible','off');
				set(0,'CurrentFigure',ih)
				plot([1 gend],[0 0],'-r','linewidth',0.5);
				hold on
				for j = 1:9
					plot(1:gend,DistribSmooth(k,:,j),'-g','linewidth',0.5)
				end
				plot(1:gend,MeanSmooth(k,:),'-k','linewidth',1)
				xlim([1 gend]);
				hold off
				if ~isempty(options_.XTick)
				  	set(gca,'XTick',options_.XTick)
				  	set(gca,'XTickLabel',options_.XTickLabel)
				end
				eval(['print -depsc2 ' M_.fname '_SmoothedVariable_' name]);
				eval(['print -dpdf ' M_.fname '_SmoothedVariable_' name]);
				saveas(ih,[M_.fname '_SmoothedVariable_' name '.fig']);
				if TeX
				  	texname = deblank(M_.endo_names_tex(oo_.dr.order_var(k),:));
				  	TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
				end    
				set(0,'CurrentFigure',hh)
				title(name,'Interpreter','none')
				if ~isempty(options_.XTick)
		  			set(gca,'XTick',options_.XTick)
		  			set(gca,'XTickLabel',options_.XTickLabel)
				end
			end
	    	eval(['print -depsc2 ' M_.fname '_SmoothedVariables' int2str(plt)]);
	    	eval(['print -dpdf ' M_.fname '_SmoothedVariables' int2str(plt)]);
	    	saveas(hh,[M_.fname '_SmoothedVariables' int2str(plt) '.fig']);
	    	if options_.nograph, close(hh), end
	    	if TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
				for jj = 1:nstar
					fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
				end    
				fprintf(fidTeX,'\\centering \n');
				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedVariables%s}\n',M_.fname,int2str(plt));
				fprintf(fidTeX,'\\caption{Smoothed variables.}');
				fprintf(fidTeX,'\\label{Fig:SmoothedVariables:%s}\n',int2str(plt));
				fprintf(fidTeX,'\\end{figure}\n');
				fprintf(fidTeX,' \n');
	    	end    
		end
	    hh = figure('Name','Smoothed variables');
	    NAMES = [];
	    TEXNAMES = [];
	    for i=1:nvar-(nbplt-1)*nstar
	      	k = (nbplt-1)*nstar+i;
	      	set(0,'CurrentFigure',hh)
	      	if lr ~= 0
				subplot(lr,lc,i);
	      	else
				subplot(nr,nc,i);
	      	end    
	      	plot([1 gend],[0 0],'-r','linewidth',0.5)
	      	hold on
	      	for j = 1:9
				plot(1:gend,DistribSmooth(k,:,j),'-g','linewidth',0.5)
	      	end
	      	plot(1:gend,MeanSmooth(k,:),'-k','linewidth',1)
	      	xlim([1 gend]);
	      	hold off
	      	name     = deblank(M_.endo_names(oo_.dr.order_var(k),:));
	      	NAMES    = strvcat(NAMES,name);
	      	ih = figure('Visible','off');
	      	set(0,'CurrentFigure',ih)
	      	plot([1 gend],[0 0],'-r','linewidth',0.5);
	      	hold on
	      	for j = 1:9
				plot(1:gend,DistribSmooth(k,:,j),'-g','linewidth',0.5)
	      	end
	      	plot(1:gend,MeanSmooth(k,:),'-k','linewidth',1)
	      	xlim([1 gend]);
	      	hold off
	      	if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
	      	end
	      	eval(['print -depsc2 ' M_.fname '_SmoothedVariable_' name]);
	      	eval(['print -dpdf ' M_.fname '_SmoothedVariable_' name]);
	      	saveas(ih,[M_.fname '_SmoothedVariable_' name '.fig']);
	      	if TeX
				texname  = deblank(M_.endo_names_tex(oo_.dr.order_var(k),:));
				TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
	      	end
	      	set(0,'CurrentFigure',hh)
	      	title(name,'Interpreter','none');
	      	if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
	      	end
		end
	    eval(['print -depsc2 ' M_.fname '_SmoothedVariables' int2str(nbplt)]);
	    eval(['print -dpdf ' M_.fname '_SmoothedVariables' int2str(nbplt)]);
	    saveas(hh,[M_.fname '_SmoothedVariables' int2str(nbplt) '.fig']);
	    if options_.nograph, close(hh), end
	    if TeX
	    	fprintf(fidTeX,'\\begin{figure}[H]\n');
	      	for jj = 1:size(NAMES,1);
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	      	end    
	      	fprintf(fidTeX,'\\centering \n');
	      	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedVariables%s}\n',M_.fname,int2str(nbplt));
	      	fprintf(fidTeX,'\\caption{Smoothed variables.}');
	      	fprintf(fidTeX,'\\label{Fig:SmoothedVariables:%s}\n',int2str(nbplt));
	      	fprintf(fidTeX,'\\end{figure}\n');
	      	fprintf(fidTeX,'\n');
	      	fprintf(fidTeX,'%% End of TeX file.\n');
	      	fclose(fidTeX);
	    end
	end % nbplt == 1 (smooth variables)	
	%%
	%% Smoothed observation error
	%%
	if nvn
		number_of_plots_to_draw = 0;
		index = [];
		for i=1:nvobs
			if max(abs(MeanError(10:end))) > 0.000000001
		    	number_of_plots_to_draw = number_of_plots_to_draw + 1;
		    	index = cat(1,index,i);
		  	end
		end
		[nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
		if TeX
		  	fidTeX = fopen([M_.fname '_SmoothedObservationErrors.TeX'],'w');
		  	fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
		  	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
		  	fprintf(fidTeX,' \n');
		end
		if nbplt == 1
		  	hh = figure('Name','Smoothed observation errors');    
		  	NAMES = [];
		  	if TeX; TEXNAMES = []; end;
		  	for i=1:number_of_plots_to_draw
		  	  	set(0,'CurrentFigure',hh)
		  	  	subplot(nr,nc,i);
		  	  	plot([1 gend],[0 0],'-r','linewidth',0.5);
		  	  	hold on
		  	  	for j = 1:9
		  	  		plot(1:gend,DistribError(:,index(i),j),'-g','linewidth',0.5)
		  	  	end
		  	  	plot(1:gend,MeanError(:,index(i)),'-k','linewidth',1)
		  	  	xlim([1 gend]);
		  	  	hold off
		  	  	ih = figure('Visible','off');
		  	  	set(0,'CurrentFigure',ih)
		  	  	plot([1 gend],[0 0],'-r','linewidth',0.5);
		  	  	hold on
		  	  	for j = 1:9
		  	  		plot(1:gend,DistribError(:,index(i),j),'-g','linewidth',0.5)
		  	  	end
		  	  	plot(1:gend,MeanError(:,index(i)),'-k','linewidth',1)
		  	  	xlim([1 gend]);
		  	  	hold off
		  	  	name    = deblank(options_.varobs(index(i),:));
		  	  	NAMES   = strvcat(NAMES,name);
		  	  	if ~isempty(options_.XTick)
		  	  	  	set(gca,'XTick',options_.XTick)
		  	  	  	set(gca,'XTickLabel',options_.XTickLabel)
		  	  	end
		  	  	eval(['print -depsc2 ' M_.fname '_SmoothedObservationError_' name]);
		  	  	eval(['print -dpdf ' M_.fname '_SmoothedObservationError_' name]);
		  	  	saveas(ih,[M_.fname '_SmoothedObservationError_' name '.fig']);
		  	  	if TeX
		  	  	  	texname = deblank(options_.varobs_TeX(index(i),:));
		  	  	  	TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
		  	  	end
		  	  	set(0,'CurrentFigure',hh)
		  	  	title(name,'Interpreter','none')
		  	  	if ~isempty(options_.XTick)
		  	  	  	set(gca,'XTick',options_.XTick)
		  	  	  	set(gca,'XTickLabel',options_.XTickLabel)
		  	  	end
		  	end
		  	eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(1)]);
		  	eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(1)]);
		  	saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(1) '.fig']);
		  	if options_.nograph, close(hh), end
		  	if TeX
		    	fprintf(fidTeX,'\\begin{figure}[H]\n');
		    	for jj = 1:M_.exo_nbr
		      		fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
		    	end    
		    	fprintf(fidTeX,'\\centering \n');
		    	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(1));
		    	fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
		    	fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(1));
		    	fprintf(fidTeX,'\\end{figure}\n');
		    	fprintf(fidTeX,' \n');
		    	fprintf(fidTeX,'%% End of TeX file.\n');
		    	fclose(fidTeX);
		  	end    
		else
			for plt = 1:nbplt-1
		    	hh = figure('Name','Smoothed observation errors');
		    	NAMES = [];
		    	if TeX; TEXNAMES = []; end;
		    	for i=1:nstar
		      		k = (plt-1)*nstar+i;
		      		set(0,'CurrentFigure',hh)
		      		subplot(nr,nc,i);
		      		plot([1 gend],[0 0],'-r','linewidth',0.5)
		      		hold on
		      		for j = 1:9
						plot(1:gend,DistribError(:,index(k),j),'-g','linewidth',0.5)
		      		end
		      		plot(1:gend,MeanError(:,index(k)),'-k','linewidth',1)
		      		xlim([1 gend]);
		      		hold off
		      		name = deblank(options_.varobs(index(k),:));
		      		NAMES = strvcat(NAMES,name);
		      		ih = figure('Visible','off');
		      		set(0,'CurrentFigure',ih)
		      		plot([1 gend],[0 0],'-r','linewidth',0.5);
		      		hold on
		      		for j = 1:9
						plot(1:gend,DistribError(:,index(k),j),'-g','linewidth',0.5)
		      		end
		      		plot(1:gend,MeanError(:,index(k)),'-k','linewidth',1)
		      		xlim([1 gend]);
		      		hold off
		      		if ~isempty(options_.XTick)
						set(gca,'XTick',options_.XTick)
						set(gca,'XTickLabel',options_.XTickLabel)
		      		end
		      		eval(['print -depsc2 ' M_.fname '_SmoothedObservationError_' name]);
		      		eval(['print -dpdf ' M_.fname '_SmoothedObservationError_' name]);
		      		saveas(ih,[M_.fname '_SmoothedObservationError_' name '.fig']);
		      		if TeX
						texname = deblank(options_.varobs_TeX(index(k),:));
						TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
		      		end    
		      		set(0,'CurrentFigure',hh)
		      		title(name,'Interpreter','none')
		      		if ~isempty(options_.XTick)
						set(gca,'XTick',options_.XTick)
						set(gca,'XTickLabel',options_.XTickLabel)
		      		end
		    	end
		    	eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(plt)]);
		    	eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(plt)]);
		    	saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(plt) '.fig']);
		    	if options_.nograph, close(hh), end
		    	if TeX
		      		fprintf(fidTeX,'\\begin{figure}[H]\n');
		      		for jj = 1:nstar
						fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
		      		end    
		      		fprintf(fidTeX,'\\centering \n');
		      		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(plt));
		      		fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
		      		fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
		      		fprintf(fidTeX,'\\end{figure}\n');
		      		fprintf(fidTeX,' \n');
		    	end    
		  	end
		  	hh = figure('Name','Smoothed observation errors');
		  	NAMES = [];
		  	if TeX; TEXNAMES = []; end;
		  	for i=1:M_.exo_nbr-(nbplt-1)*nstar
		    	k = (nbplt-1)*nstar+i;
		    	set(0,'CurrentFigure',hh)
		    	if lr ~= 0
		    		subplot(lr,lc,i);
		    	else
		    		subplot(nr,nc,i);
		    	end    
		    	plot([1 gend],[0 0],'-r','linewidth',0.5)
		    	hold on
		    	for j = 1:9
		     		plot(1:gend,DistribError(:,index(k),j),'-g','linewidth',0.5)
		    	end
		    	plot(1:gend,MeanError(:,index(k)),'-k','linewidth',1)
		    	xlim([1 gend]);
		    	hold off
		    	name = deblank(options_.varobs(index(k),:));
		    	NAMES = strvcat(NAMES,name);
		    	ih = figure('Visible','off');
		    	set(0,'CurrentFigure',ih)
		    	plot([1 gend],[0 0],'-r','linewidth',0.5);
		    	hold on
		    	for j = 1:9
		      		plot(1:gend,DistribError(:,index(k),j),'-g','linewidth',0.5)
		    	end
		    	plot(1:gend,MeanError(:,index(k)),'-k','linewidth',1)
		    	xlim([1 gend]);
		    	hold off
		    	if ~isempty(options_.XTick)
		      		set(gca,'XTick',options_.XTick)
		      		set(gca,'XTickLabel',options_.XTickLabel)
		    	end
		    	eval(['print -depsc2 ' M_.fname '_SmoothedObservationError_' name]);
		    	eval(['print -dpdf ' M_.fname '_SmoothedObservationError_' name]);
		    	saveas(ih,[M_.fname '_SmoothedObservationError_' name '.fig']);
		    	if TeX
		      		texname  = deblank(options_.varobs_TeX(index(k),:));
		      		TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
		    	end
		    	set(0,'CurrentFigure',hh)
		    	title(name,'Interpreter','none');
		    	if ~isempty(options_.XTick)
		      		set(gca,'XTick',options_.XTick)
		      		set(gca,'XTickLabel',options_.XTickLabel)
		    	end
		  	end
		  	eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(nbplt)]);
		  	eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(nbplt)]);
		  	saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(nbplt) '.fig']);
		  	if options_.nograph, close(hh), end
		  	if TeX
		    	fprintf(fidTeX,'\\begin{figure}[H]\n');
		    	for jj = 1:size(NAMES,1);
		      		fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
		    	end    
		    	fprintf(fidTeX,'\\centering \n');
		    	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(nbplt));
		    	fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
		    	fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(nbplt));
		    	fprintf(fidTeX,'\\end{figure}\n');
		    	fprintf(fidTeX,'\n');
		    	fprintf(fidTeX,'%% End of TeX file.\n');
		    	fclose(fidTeX);
		  	end
		end % nbplt == 1 (smooth observation errors)
	end		
	%%
	%%	Historical and smoothed variabes
	%%
	[atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(post_mean,gend,data);
	yf = zeros(gend,nvobs,B);
	if options_.prefilter == 1
		yf = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend)+trend_coeff*[0:gend-1];
	elseif options_.loglinear == 1
		yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
		     trend_coeff*[0:gend-1];
	else
		yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
		     trend_coeff*[0:gend-1];
	end
	[nbplt,nr,nc,lr,lc,nstar] = pltorg(nvobs);
	if TeX
	  	fidTeX = fopen([M_.fname '_HistoricalAndSmoothedVariables.TeX'],'w');
	  	fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
	  	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
	  	fprintf(fidTeX,' \n');
	end    
	if nbplt == 1
		hh = figure('Name','Historical and smoothed variables');
		NAMES = [];
		if TeX ; TEXNAMES = []; end;
		for i=1:nvobs
			subplot(nr,nc,i);
			plot(1:gend,yf(i,1:end),'-r','linewidth',1) %yf(i,2:end)
			hold on
			plot(1:gend,rawdata(:,i),'-k','linewidth',1)
			xlim([1 gend]);
			hold off
			name    = options_.varobs(i,:);
			NAMES   = strvcat(NAMES,name);
			if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
			end
			if TeX
				texname = options_.varobs_TeX(i,1);
			  	TEXNAMES   = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
			end
			title(name,'Interpreter','none')
		end
		eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1)]);
		eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1)]);
		saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(1) '.fig']);
		if options_.nograph, close(hh), end
		if TeX
			fprintf(fidTeX,'\\begin{figure}[H]\n');
			for jj = 1:n_varobs
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
			end    	
			fprintf(fidTeX,'\\centering \n');
			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(1));
			fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
			fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(1));
			fprintf(fidTeX,'\\end{figure}\n');
			fprintf(fidTeX,'\n');
			fprintf(fidTeX,'%% End of TeX file.\n');
			fclose(fidTeX);
		end    
	else
		for plt = 1:nbplt-1
			hh = figure('Name','Historical and smoothed variables');
			set(0,'CurrentFigure',hh)
			NAMES = [];
			if TeX; TEXNAMES = []; end;
			for i=1:nstar
				k = (plt-1)*nstar+i;
			  	subplot(nr,nc,i);
			  	plot(1:gend,yf(k,2:end),'-r','linewidth',1)
			  	hold on
			  	plot(1:gend,rawdata(:,k),'-k','linewidth',1)
			  	xlim([1 gend]);
			  	hold off
			  	name = options_.varobs(k,:);
			  	NAMES = strvcat(NAMES,name);
			  	if ~isempty(options_.XTick)
			    	set(gca,'XTick',options_.XTick)
			    	set(gca,'XTickLabel',options_.XTickLabel)
			  	end
			  	if TeX
			    	texname = options_.varobs_TeX(k,:);
			    	TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
			  	end    
			  	title(name,'Interpreter','none')
			end
			eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)]);
			eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)]);
			saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(plt) '.fig']);
			if options_.nograph, close(hh), end
			if TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
			  	for jj = 1:nstar
			    	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
			  	end    	
			  	fprintf(fidTeX,'\\centering \n');
			  	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(plt));
			  	fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
			  	fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
			  	fprintf(fidTeX,'\\end{figure}\n');
			  	fprintf(fidTeX,'\n');
			end    
		end
		hh = figure('Name','Historical and smoothed variables');
		set(0,'CurrentFigure',hh)
		NAMES = [];
		if TeX; TEXNAMES = []; end;
		for i=1:nobs-(nbplt-1)*nstar
			k = (nbplt-1)*nstar+i;
			if lr ~= 0
				subplot(lr,lc,i);
			else
				subplot(nr,nc,i);
			end    
			plot(1:gend,yf(k,2:end),'-r','linewidth',1)
			hold on
			plot(1:gend,rawdata(:,k),'-k','linewidth',1)
			xlim([1 gend]);
			hold off
			name = options_.varobs(k,:);
			NAMES    = strvcat(NAMES,name);
			if ~isempty(options_.XTick)
			  	set(gca,'XTick',options_.XTick)
			  	set(gca,'XTickLabel',options_.XTickLabel)
			end
			if TeX
			  	texname  = options_.varobs_TeX(k,:);
			  	TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(texname) ' $']);
			end
			title(name,'Interpreter','none');
		end
		eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
		eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
		saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt) '.fig']);
		if options_.nograph, close(hh), end
		if TeX
			fprintf(fidTeX,'\\begin{figure}[H]\n');
			for jj = 1:size(NAMES,1)
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
			end    	
			fprintf(fidTeX,'\\centering \n');
			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(nbplt));
			fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
			fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(nbplt));
			fprintf(fidTeX,'\\end{figure}\n');
			fprintf(fidTeX,'\n');
			fprintf(fidTeX,'%% End of TeX file.\n');
			fclose(fidTeX);
		end    
	end
end % options_.smoother
if options_.filtered_vars
	fprintf('MH: Filtered variables...\n')
	MeanFilter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	MedianFilter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	StdFilter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend);
	DistribFilter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,9);
	HPDFilter = zeros(size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic,gend,2);
	for i = 1:size(oo_.dr.ghx,2)+oo_.dr.nfwrd+oo_.dr.nstatic;
		for t = 1:gend
			StartLine = 0;
		  	for file = 1:sfil_filt;
		    	instr = [M_.fname '_filter' int2str(file)];
		    	eval(['load ' instr]);
		    	MeanFilter(i,t) = MeanFilter(i,t)+sum(stock_filter(i,t,:),3);
		    	DeProfundis = size(stock_filter,3); 
		    	tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_filter(i,t,:)); 
		    	StartLine = StartLine+DeProfundis;
		  	end
		  	tmp = sort(tmp);
		  	MedianFilter(i,t) = tmp(round(B*0.5));
		  	StdFilter(i,t) = std(tmp);
		  	DistribFilter(i,t,:) = reshape(tmp(deciles),1,1,9);
		  	tt = floor(options_.mh_conf_sig*B);
		  	a = 1; 
		  	b = tt;
		  	tmp2 = [1;tt;tmp(tt)-tmp(1)];
		  	while b <= B
		    	tmp1 = [a;b;tmp(b)-tmp(a)];
		    	a = a + 1;
		    	b = b + 1;
		    	if tmp1(3,1) < tmp2(3,1)
		      		tmp2 = tmp1;     
		    	end    
		  	end
		  	HPDFilter(i,t,1) = tmp(tmp2(1,1));
		  	HPDFilter(i,t,2) = tmp(tmp2(2,1));
		end
		disp(['    Variable: ' deblank(M_.endo_names(oo_.dr.order_var(i),:))]);	
	end
    clear stock_filter;
    MeanFilter = MeanFilter/B;
    for i=1:size(M_.endo_names,1)
		eval(['oo_.PosteriorFilteredVariables.Mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = MeanFilter(i,:)'';']);
		eval(['oo_.PosteriorFilteredVariables.Median.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = MedianFilter(i,:)'';']);
		eval(['oo_.PosteriorFilteredVariables.Std.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = StdFilter(i,:)'';']);
		eval(['oo_.PosteriorFilteredVariables.Distribution.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(DistribFilter(i,:,:))'';']);
		eval(['oo_.PosteriorFilteredVariables.HPDinf.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(HPDFilter(i,:,1))'';']);
		eval(['oo_.PosteriorFilteredVariables.HPDsup.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = squeeze(HPDFilter(i,:,2))'';']);
	end
    fprintf('MH: Filtered variables, done!\n')
    disp(' ')
end	       
%%
%% 	Posterior IRFs. Instead of displaying the IRFs associated to the posterior mean
%%	of the structural parameters (by calling stoch_simul after estimation), 
%%	metropolis.m will display the posterior mean of the IRFs and the deciles of 
%%	the IRFs' posterior distribution. All the results are saved in the global 
%%	structure oo_ (posterior medians, posterior standard deviations and posterior HPD   
%%	intervals are also computed and saved).
%%
if options_.irf
	if B <= MAX_nirfs
		stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,B);
	elseif nvn & B > MAX_nirfs
	  	stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,MAX_nirfs);
	end
	h = waitbar(0,'Bayesian IRFs...');
	if nfile-ffil+1>1
	   	sfil_irf = 1;
	   	irun_irf = 0;
	   	for b = 1:B;
	     	irun_irf = irun_irf+1;
	     	tmp = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr);
	     	choose_an_mh_file = rand;
	     	mh_file_number = ...
		 	FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	     	if isempty(mh_file_number)
	     		mh_file_number = ffil;
	     	else    
	     	 	mh_file_number = mh_file_number(1);
	     	end    
	     	eval(['load ' instr1 int2str(mh_file_number) instr2]);
	     	clear post2 logpo2;
	     	deep  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	     	offset = nvx+nvn+ncx+ncn;
	     	for i=1:estim_params_.np
	     		assignin('base',deblank(estim_params_.param_names(i,:)),deep(i+offset));
	     	end
	     	dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
	     	if nvx
	     		ip = 1;
	     	  	for i=1:nvx
		 			k = estim_params_.var_exo(i,1);
		 			M_.Sigma_e(k,k) = deep(ip)*deep(ip);
		 			ip = ip+1;
	     	  	end
	     	end
	     	if ncx
	     		ip = nvx+nvn+1;
	     	  	for i=1:ncx
		 			k1 = estim_params_.corrx(i,1);
		 			k2 = estim_params_.corrx(i,2);
		 			M_.Sigma_e(k1,k2) = deep(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
		 			M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
		 			ip = ip+1;
	     	  	end
	     	end
	     	SS(M_.exo_name_orig_ord,M_.exo_name_orig_ord)=M_.Sigma_e+1e-14* ...
		    eye(M_.exo_nbr);
	     	SS = transpose(chol(SS));
	     	tit(M_.exo_name_orig_ord,:) = M_.exo_names;
	     	for i = 1:M_.exo_nbr
	     		if SS(i,i) > 1e-13
		 			y=irf(dr,SS(M_.exo_name_orig_ord,i), options_.irf, options_.drop, ...
		 	      		options_.replic, options_.order);
					if options_.relative_irf
					  y = 100*y/cs(i,i); 
					end

		 			for j = 1:size(M_.endo_names,1)
		 	  			if max(y(j,:)) - min(y(j,:)) > 1e-10 
		 	    			tmp(:,j,i) = transpose(y(j,:));
		 	  			end	
		 			end	
	     	  	end
	     	end
	     	if irun_irf < MAX_nirfs
	       		stock_irf(:,:,:,irun_irf) = tmp;
	     	else
	       		stock_irf(:,:,:,irun_irf) = tmp;
	       		instr = [M_.fname '_irf' int2str(sfil_irf) ' stock_irf;'];
	       		eval(['save ' instr]);
	       		sfil_irf = sfil_irf + 1;
	       		irun_irf = 0;
	       		stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,MAX_nirfs);
	     	end	
	     	waitbar(b/B,h);    
	   	end
	   	clear tmp;
	   	if irun_irf
	    	stock_irf = stock_irf(:,:,:,1:irun_irf);
	     	instr = [M_.fname '_irf' int2str(sfil_irf) ' stock_irf;'];
	     	eval(['save ' instr]);
	   	end
	   	clear stock_irf;
	else		
		sfil_irf = 1;
	   	irun_irf = 0;
	   	eval(['load ' instr1 int2str(ffil) instr2]);
	   	NumberOfSimulations = length(logpo2);
	   	clear post2 logpo2;
	   	for b = 1:B;
	    	irun_irf = irun_irf+1;
	     	tmp = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr);
	     	deep  = x2(floor(rand*NumberOfSimulations)+1,:);
	     	offset = nvx+nvn+ncx+ncn;
	     	for i=1:estim_params_.np
	       		assignin('base',deblank(estim_params_.param_names(i,:)),deep(i+offset));
	     	end
	     	dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
	     	if nvx
	       		ip = 1;
	       		for i=1:nvx
		 			k = estim_params_.var_exo(i,1);
		 			M_.Sigma_e(k,k) = deep(ip)*deep(ip);
		 			ip = ip+1;
	       		end
	     	end
	     	if ncx
	       		ip = nvx+nvn+1;
	       		for i=1:ncx
		 			k1 = estim_params_.corrx(i,1);
		 			k2 = estim_params_.corrx(i,2);
		 			M_.Sigma_e(k1,k2) = deep(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
		 			M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
		 			ip = ip+1;
	       		end
	     	end
	     	SS(M_.exo_name_orig_ord,M_.exo_name_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
	     	SS = transpose(chol(SS));
	     	tit(M_.exo_name_orig_ord,:) = M_.exo_names;
	     	for i = 1:M_.exo_nbr
	       		if SS(i,i) > 1e-13
		 			y=irf(oo_.dr,SS(M_.exo_name_orig_ord,i), options_.irf, options_.drop, ...
		       			options_.replic, options_.order);
					if options_.relative_irf
					  y = 100*y/cs(i,i); 
					end
		 			for j = 1:size(M_.endo_names,1)
		   				if max(y(j,:)) - min(y(j,:)) > 1e-10 
		     				tmp(:,j,i) = transpose(y(j,:));
		   				end	
		 			end	
	       		end
			end
	    	if irun_irf < MAX_nirfs
	    		stock_irf(:,:,:,irun_irf) = tmp;
	    	else
	    		stock_irf(:,:,:,irun_irf) = tmp;
	       		instr = [M_.fname '_irf' int2str(sfil_irf) ' stock_irf;'];
	       		eval(['save ' instr]);
	       		sfil_irf = sfil_irf + 1;
	       		irun_irf = 0;
	       		stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,MAX_nirfs);
			end	
	    	waitbar(b/B,h);    
		end
	   	if irun_irf
	    	stock_irf = stock_irf(:,:,:,1:irun_irf);
	     	instr = [M_.fname '_irf' int2str(sfil_irf) ' stock_irf;'];
	     	eval(['save ' instr]);
	   	end
	   	clear stock_irf;
	end		
	close(h)
	%%
	%% 	Now i compute some statistics (mean, median, std, deciles, HPD intervals)
	%%
	tmp = zeros(B,1);
	MeanIRF = zeros(options_.irf,nvar,M_.exo_nbr);
	MedianIRF = zeros(options_.irf,nvar,M_.exo_nbr);
	StdIRF = zeros(options_.irf,nvar,M_.exo_nbr);
	DistribIRF = zeros(options_.irf,nvar,M_.exo_nbr,9);
	HPDIRF = zeros(options_.irf,nvar,M_.exo_nbr,2);
	fprintf('MH: Posterior IRFs...\n')
	for i = 1:M_.exo_nbr
		for j = 1:nvar
	    	for k = 1:options_.irf
	       		StartLine = 0;
	       		for file = 1:sfil_irf;
		 			instr = [M_.fname '_irf' int2str(file)];
		 			eval(['load ' instr]);
		 			MeanIRF(k,j,i) = MeanIRF(k,j,i)+sum(stock_irf(k,SelecVariables(j),i,:),4);
		 			DeProfundis = size(stock_irf,4); 
		 			tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_irf(k,SelecVariables(j),i,:)); 
		 			StartLine = StartLine+DeProfundis;
	       		end
	       		tmp = sort(tmp);
	       		MedianIRF(k,j,i) = tmp(round(B*0.5));
	       		StdIRF(k,j,i) = std(tmp);
	       		DistribIRF(k,j,i,:) = reshape(tmp(deciles),1,1,1,9);
	       		tt = floor(options_.mh_conf_sig*B);
	       		a = 1; 
	       		b = tt;
	       		tmp2 = [1;tt;tmp(tt)-tmp(1)];
	       		while b <= B
		 			tmp1 = [a;b;tmp(b)-tmp(a)];
		 			a = a + 1;
		 			b = b + 1;
		 			if tmp1(3,1) < tmp2(3,1)
		   				tmp2 = tmp1;     
		 			end    
	       		end
	       		HPDIRF(k,j,i,1) = tmp(tmp2(1,1));
	       		HPDIRF(k,j,i,2) = tmp(tmp2(2,1));
	     	end
	     	disp(['    Variable: ' deblank(M_.endo_names(SelecVariables(j),:)) ', orthogonalized shock to ' deblank(tit(i,:))])	
	   	end	
	end
	clear stock_irf;
	MeanIRF = MeanIRF/B;
	for i = 1:M_.exo_nbr
		for j = 1:nvar
	    	eval(['oo_.PosteriorIRF.Mean.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = MeanIRF(:,j,i);']);
	     	eval(['oo_.PosteriorIRF.Median.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = MedianIRF(:,j,i);']);
	     	eval(['oo_.PosteriorIRF.Std.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = StdIRF(:,j,i);']);
	     	eval(['oo_.PosteriorIRF.Distribution.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(DistribIRF(:,j,i,:));']);
	     	eval(['oo_.PosteriorIRF.HPDinf.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF(:,j,i,1));']);
			eval(['oo_.PosteriorIRF.HPDsup.' deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF(:,j,i,2));']);
	   	end
	end	
	%%
	%% 	Finally i build the plots.
	%%
	if TeX
		fidTeX = fopen([M_.fname '_BayesianIRF.TeX'],'w');
		fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
		fprintf(fidTeX,' \n');
	end
	tit(M_.exo_name_orig_ord,:) = M_.exo_names;
	if TeX; titTeX(M_.exo_name_orig_ord,:) = M_.exo_names_tex; end;
	for i=1:M_.exo_nbr
		number_of_plots_to_draw = 0;
	   	index = [];
	   	for j=1:nvar
	    	if MeanIRF(1,j,i)
	       		number_of_plots_to_draw = number_of_plots_to_draw + 1;
	       		index = cat(1,index,j);
	     	end
	   	end
	   	[nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);	
	   	if nbplt == 1
		  if options_.relative_irf
		    hh = figure('Name',['Relative response to orthogonalized' ...
					' shock to ' tit(i,:)]);
		  else
		    hh = figure('Name',['Orthogonalized shock to ' tit(i, ...
						  :)]);
		  end
	     	NAMES = [];
	     	if TeX; TEXNAMES = []; end;
	     	for j=1:number_of_plots_to_draw
	       		set(0,'CurrentFigure',hh)
	       		subplot(nr,nc,j);
	       		plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
	       		hold on
	       		for k = 1:9
		 			plot(1:options_.irf,DistribIRF(:,index(j),i,k),'-g','linewidth',0.5)
	       		end
	       		plot(1:options_.irf,MeanIRF(:,index(j),i),'-k','linewidth',1)
	       		xlim([1 options_.irf]);
	       		hold off
	       		name    = deblank(M_.endo_names(SelecVariables(index(j)),:));
	       		NAMES   = strvcat(NAMES,name);
	       		if TeX
		 			texname = deblank(M_.endo_names_tex(SelecVariables(index(j)),:));
		 			TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
	       		end
	       		title(name,'Interpreter','none')
	     	end
	     	eval(['print -depsc2 ' M_.fname '_Bayesian_IRF_' deblank(tit(i,:))]);
	     	eval(['print -dpdf ' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:))]);
	     	saveas(hh,[M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) '.fig']);
	     	if options_.nograph, close(hh), end
	     	if TeX
	       		fprintf(fidTeX,'\\begin{figure}[H]\n');
	       		for jj = 1:number_of_plots_to_draw
		 			fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	       		end    
	       		fprintf(fidTeX,'\\centering \n');
	       		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRF_%s}\n',M_.fname,deblank(tit(i,:)));
			if options_.relative_irf
			  fprintf(fidTeX,['\\caption{Bayesian relative' ...
					  ' IRF.}']);
			else
			    fprintf(fidTeX,'\\caption{Bayesian IRF.}');
			end
	       		fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s}\n',deblank(tit(i,:)));
	       		fprintf(fidTeX,'\\end{figure}\n');
	       		fprintf(fidTeX,' \n');
	    	end    
		elseif nbplt > 1
	    	for fig = 1:nbplt-1
		  if options_.relative_irf
		    hh = figure('Name',['Relative response to orthogonalized' ...
					' shock to ' tit(i,:) ' figure ' int2str(fig) '.']);
		  else
		    hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ...
					' figure ' int2str(fig) '.']);
		  end
		  NAMES = [];
	       		if TeX; TEXNAMES = []; end;
	       		for j=1:nstar
		 			jj = (fig-1)*nstar + j;
		 			subplot(nr,nc,j);
		 			plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
		 			hold on
		 			for k = 1:9
		   				plot(1:options_.irf,DistribIRF(:,index(jj),i,k),'-g','linewidth',0.5)
		 			end
		 			plot(1:options_.irf,MeanIRF(:,index(jj),i),'-k','linewidth',1)
		 			xlim([1 options_.irf]);
		 			hold off
		 			name    = deblank(M_.endo_names(SelecVariables(index(jj)),:));
		 			NAMES   = strvcat(NAMES,name);
		 			if TeX
		   				texname = deblank(M_.endo_names_tex(SelecVariables(index(jj)),:));
		   				TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
		 			end
		 			title(name,'Interpreter','none')
	       		end
	       		eval(['print -depsc2 ' M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) int2str(fig)]);
	       		eval(['print -dpdf ' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) int2str(fig)]);
	       		saveas(hh,[M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) int2str(fig) '.fig']);
	       		if options_.nograph, close(hh), end
	       		if TeX
			  fprintf(fidTeX,'\\begin{figure}[H]\n');
			  for jj = 1:nstar
			    fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
			  end    
			  fprintf(fidTeX,'\\centering \n');
			  fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRF_%s%s}\n',M_.fname,deblank(tit(i,:)),int2str(fig));
			  if options_.relative_irf == 1
			    fprintf(fidTeX,['\\caption{Bayesian relative' ...
					    ' IRF.}']);
			  else
			    fprintf(fidTeX,'\\caption{Bayesian IRF.}');
			  end
			  fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s:%s}\n',deblank(tit(i,:)), int2str(fig));
			  fprintf(fidTeX,'\\end{figure}\n');
			  fprintf(fidTeX,' \n');
			end    
		end
		hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ' figure ' int2str(nbplt) '.']);
		NAMES = [];
		if TeX; TEXNAMES = []; end;
		for j=1:number_of_plots_to_draw -(nbplt-1)*nstar
		  jj = (nbplt-1)*nstar + j;
		  subplot(nr,nc,j);
		  plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
		  hold on
		  for k = 1:9
		    plot(1:options_.irf,DistribIRF(:,index(jj),i,k),'-g','linewidth',0.5)
		  end
		  plot(1:options_.irf,MeanIRF(:,index(jj),i),'-k','linewidth',1)
		  xlim([1 options_.irf]);
		  hold off
		  name    = deblank(M_.endo_names(SelecVariables(index(jj)),:));
		  NAMES   = strvcat(NAMES,name);
		  if TeX
		    texname = deblank(M_.endo_names_tex(SelecVariables(index(jj)),:));
		    TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
		  end
		  title(name,'Interpreter','none')
		end
		eval(['print -depsc2 ' M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) int2str(nbplt)]);
		eval(['print -dpdf ' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) int2str(nbplt)]);
		saveas(hh,[M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) int2str(nbplt) '.fig']);
		if options_.nograph, close(hh), end
		if TeX
		  fprintf(fidTeX,'\\begin{figure}[H]\n');
		  for jj = 1:nstar
		    fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
		  end    
		  fprintf(fidTeX,'\\centering \n');
		  fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRF_%s%s}\n',M_.fname,deblank(tit(i,:)),int2str(nbplt));
		  fprintf(fidTeX,'\\caption{Bayesian IRF.}');
		  fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s:%s}\n',deblank(tit(i,:)), int2str(nbplt));
		  fprintf(fidTeX,'\\end{figure}\n');
		  fprintf(fidTeX,' \n');
		end    
		else % nbplt = 0
		  disp('There''s nothing to plot here!')
		end
	end
	if TeX
	  fprintf(fidTeX,'%% End of TeX file.\n');
	  fclose(fidTeX);
	end
	fprintf('MH: Posterior IRFs, done!\n');
end
%%
%% 	Posterior theoretical moments. Instead of displaying the posterior moments 
%%  associated to the posterior mean of the structural parameters (by calling 
%%  stoch_simul after estimation), metropolis.m will display the posterior mean 
%%  of the theoretical moments and the posterior HPD intervals of theoretical
%%  moments. All the results are saved in the global structure oo_ (posterior 
%%  medians, posterior standard deviations and posterior deciles are also
%%	computed and saved).
%%
if options_.moments_varendo
	if ~isempty(options_.unit_root_vars)
    	vartan = []; 
		for i=1:nvar
			if isempty(strmatch(deblank(varlist(i,:)),options_.unit_root_vars,'exact'))		
				vartan = strvcat(vartan,varlist(i,:));
			end	
		end
    	varlist = vartan;
    	nvar	= size(varlist,1);
    	ivar = zeros(nvar,1);
    	for i = 1:nvar
			ivar(i) = strmatch(varlist(i,:),M_.endo_names,'exact');
		end
	else
    	nvar	= size(varlist,1);
		ivar = zeros(nvar,1);
		for i = 1:nvar
			ivar(i) = strmatch(varlist(i,:),M_.endo_names,'exact');
		end
	end
	nar = options_.ar;
	if B <= MAX_nthm1
    	stock_thm1 = zeros(nvar,B);
	elseif B > MAX_nthm1
    	stock_thm1 = zeros(nvar,MAX_nthm1);
	end
	if B <= MAX_nthm2
    	stock_thm2 = zeros(nvar,nvar,B);
	elseif B > MAX_nthm2
    	stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
	end
	if B <= MAX_nthm3
    	stock_thm3 = zeros(nvar,M_.exo_nbr,B);
	elseif B > MAX_nthm3
    	stock_thm3 = zeros(nvar,M_.exo_nbr,MAX_nthm3);
	end
	if B <= MAX_nthm4
    	stock_thm4 = zeros(nvar,nar,B);
	elseif B > MAX_nthm4
    	stock_thm4 = zeros(nvar,nar,MAX_nthm4);
	end
	h = waitbar(0,'Posterior theoretical moments...');
	if nfile-ffil+1>1
    	sfil_thm1 = 1;
		irun_thm1 = 0;
		sfil_thm2 = 1;
		irun_thm2 = 0;
    	sfil_thm3 = 1;
    	irun_thm3 = 0;
    	sfil_thm4 = 1;
    	irun_thm4 = 0;
    	for b = 1:B;
      		irun_thm1 = irun_thm1+1;
      		irun_thm2 = irun_thm2+1;
      		irun_thm3 = irun_thm3+1;
      		irun_thm4 = irun_thm4+1;
      		choose_an_mh_file = rand;
      		mh_file_number = ...
	  		FLN(find(choose_an_mh_file>=FLN(:,3)),1);
			if isempty(mh_file_number)
				mh_file_number = ffil;
			else
				mh_file_number = mh_file_number(1);
			end
			eval(['load ' instr1 int2str(mh_file_number) instr2]);
			clear post2 logpo2;
			deep  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
			offset = nvx+nvn+ncx+ncn;
			for i=1:estim_params_.np
				assignin('base',deblank(estim_params_.param_names(i,:)),deep(i+offset));
			end
			dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
			if nvx
				ip = 1;
				for i=1:nvx
	  				k = estim_params_.var_exo(i,1);
	  				M_.Sigma_e(k,k) = deep(ip)*deep(ip);
	  				ip = ip+1;
				end
			end
			if ncx
				ip = nvx+nvn+1;
				for i=1:ncx
	  				k1 = estim_params_.corrx(i,1);
	  				k2 = estim_params_.corrx(i,2);
	  				M_.Sigma_e(k1,k2) = deep(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
	  				M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
	  				ip = ip+1;
				end
			end
			if nvn
				Sigma_m = zeros(size(options_.varobs,1));
				ip = nvx+1;
				for i=1:nvn
					k = estim_params_.var_endo(i,1);
					Sigma_m(k,k) = deep(ip)*deep(ip);
					ip = ip + 1;
				end
			end
			if ncn
				ip = nvx+nvn+ncx+1;
				for i=1:ncn
	  				k1 = estim_params_.corrn(i,1);
	  				k2 = estim_params_.corrn(i,2);
	  				Sigma_m(k1,k2) = deep(ip)*sqrt(Sigma_m(k1,k1)*Sigma_m(k2,k2));
	  				Sigma_m(k2,k1) = Sigma_m(k1,k2);
	  				ip = ip+1;
				end
			end
			Gamma_y = th_autocovariances(dr,ivar);
			if options_.order == 2
				m_mean = dr.ys(ivar) + Gamma_y{options_.ar+3};
			else
				m_mean = dr.ys(ivar);
			end
			variance =  Gamma_y{1};
			if irun_thm1 < MAX_nthm1
				stock_thm1(:,irun_thm1) = m_mean;
			else
				stock_thm1(:,irun_thm1) = m_mean;
				instr = [M_.fname '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
				eval(['save ' instr]);
				sfil_thm1 = sfil_thm1 + 1;
				irun_thm1 = 0;
				stock_thm1 = zeros(nvar,MAX_nthm1);
			end
			if irun_thm2 < MAX_nthm2
				stock_thm2(:,:,irun_thm2) = variance;
			else
				stock_thm2(:,:,irun_thm2) = variance;
				instr = [M_.fname '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
				eval(['save ' instr]);
				sfil_thm2 = sfil_thm2 + 1;
				irun_thm2 = 0;
				stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
			end
			if irun_thm3 < MAX_nthm3
				stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
			else
				stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
				instr = [M_.fname '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
				eval(['save ' instr]);
				sfil_thm3 = sfil_thm3 + 1;
				irun_thm3 = 0;
				stock_thm3 = zeros(nvar,M_.exo_nbr,MAX_nthm3);
			end
			if irun_thm4 < MAX_nthm4
				for lag = 1:nar
					stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
				end	
			else
				for lag = 1:nar
					stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
				end	
				instr = [M_.fname '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
				eval(['save ' instr]);
				sfil_thm4 = sfil_thm4 + 1;
				irun_thm4 = 0;
				stock_thm4 = zeros(nvar,nar,MAX_nthm4);
			end
			waitbar(b/B,h);    
		end
		clear m_mean variance Gamma_y;
		if irun_thm1
      		stock_thm1 = stock_thm1(:,1:irun_thm1);
      		instr = [M_.fname '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
      		eval(['save ' instr]);
		end
		clear stock_thm1;
		if irun_thm2
      		stock_thm2 = stock_thm2(:,:,1:irun_thm2);
      		instr = [M_.fname '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
      		eval(['save ' instr]);
		end
		clear stock_thm2;
		if irun_thm3
			stock_thm3 = stock_thm3(:,:,1:irun_thm3);
      		instr = [M_.fname '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
      		eval(['save ' instr]);
		end
		clear stock_thm3;
		if irun_thm4
      		stock_thm4 = stock_thm4(:,:,1:irun_thm4);
      		instr = [M_.fname '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
      		eval(['save ' instr]);
		end
		clear stock_thm4;
	else		
    	sfil_thm1 = 1;
    	irun_thm1 = 0;
    	sfil_thm2 = 1;
    	irun_thm2 = 0;
    	sfil_thm3 = 1;
    	irun_thm3 = 0;
    	sfil_thm4 = 1;
    	irun_thm4 = 0;
    	eval(['load ' instr1 int2str(ffil) instr2]);
    	NumberOfSimulations = length(logpo2);
    	clear post2 logpo2;
    	for b = 1:B;
      		irun_thm1 = irun_thm1+1;
      		irun_thm2 = irun_thm2+1;
      		irun_thm3 = irun_thm3+1;
      		irun_thm4 = irun_thm4+1;
      		deep  = x2(floor(rand*NumberOfSimulations)+1,:);
      		offset = nvx+nvn+ncx+ncn;
      		for i=1:estim_params_.np
				assignin('base',deblank(estim_params_.param_names(i,:)),deep(i+offset));
			end
			dr = resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);
			if nvx
				ip = 1;
				for i=1:nvx
	  				k = estim_params_.var_exo(i,1);
	  				M_.Sigma_e(k,k) = deep(ip)*deep(ip);
	  				ip = ip+1;
				end
			end
			if ncx
				ip = nvx+nvn+1;
				for i=1:ncx
					k1 = estim_params_.corrx(i,1);
	  				k2 = estim_params_.corrx(i,2);
	  				M_.Sigma_e(k1,k2) = deep(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
	  				M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
	  				ip = ip+1;
				end
			end
			if nvn
				Sigma_m = zeros(size(options_.varobs,1));
				ip = nvx+1;
				for i=1:nvn
	  				k = estim_params_.var_endo(i,1);
	  				Sigma_m(k,k) = deep(ip)*deep(ip);
	  				ip = ip + 1;
				end
			end
			if ncn
				ip = nvx+nvn+ncx+1;
				for i=1:ncn
	  				k1 = estim_params_.corrn(i,1);
	  				k2 = estim_params_.corrn(i,2);
	  				Sigma_m(k1,k2) = deep(ip)*sqrt(Sigma_m(k1,k1)*Sigma_m(k2,k2));
	  				Sigma_m(k2,k1) = Sigma_m(k1,k2);
	  				ip = ip+1;
				end
			end
			Gamma_y = th_autocovariances(dr,ivar);
			if options_.order == 2
				m_mean = dr.ys(ivar) + Gamma_y{options_.ar+3};
			else
				m_mean = dr.ys(ivar);
			end
			variance = Gamma_y{1};
			if irun_thm1 < MAX_nthm1
				stock_thm1(:,irun_thm1) = m_mean;
			else
				stock_thm1(:,irun_thm1) = m_mean;
				instr = [M_.fname '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
				eval(['save ' instr]);
				sfil_thm1 = sfil_thm1 + 1;
				irun_thm1 = 0;
				stock_thm1 = zeros(nvar,MAX_nthm1);
			end
			if irun_thm2 < MAX_nthm2
				stock_thm2(:,:,irun_thm2) = variance;
			else
				stock_thm2(:,:,irun_thm2) = variance;
				instr = [M_.fname '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
				eval(['save ' instr]);
				sfil_thm2 = sfil_thm2 + 1;
				irun_thm2 = 0;
				stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
			end
			if irun_thm3 < MAX_nthm3
				stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
			else
				stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
				instr = [M_.fname '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
				eval(['save ' instr]);
				sfil_thm3 = sfil_thm3 + 1;
				irun_thm3 = 0;
				stock_thm3 = zeros(nvar,M_.exo_nbr,MAX_nthm3);
			end
			if irun_thm4 < MAX_nthm4
				for lag = 1:nar
					stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
				end	
			else
				for lag = 1:nar
					stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
				end
				instr = [M_.fname '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
				eval(['save ' instr]);
				sfil_thm4 = sfil_thm4 + 1;
				irun_thm4 = 0;
				stock_thm4 = zeros(nvar,nar,MAX_nthm4);
			end
			waitbar(b/B,h);    
		end
		clear m_mean variance Gamma_y;
    	if irun_thm1
    	  	stock_thm1 = stock_thm1(:,1:irun_thm1);
    	  	instr = [M_.fname '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
    	  	eval(['save ' instr]);
		end
    	clear stock_thm1;
    	if irun_thm2
    	  	stock_thm2 = stock_thm2(:,:,1:irun_thm2);
    	  	instr = [M_.fname '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
    	  	eval(['save ' instr]);
		end
    	clear stock_thm2;
    	if irun_thm3
    	  	stock_thm3 = stock_thm3(:,:,1:irun_thm3);
    	  	instr = [M_.fname '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
    	  	eval(['save ' instr]);
		end
    	clear stock_thm3;
    	if irun_thm4
    	  	stock_thm4 = stock_thm4(:,:,1:irun_thm4);
    	  	instr = [M_.fname '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
    	  	eval(['save ' instr]);
		end
    	clear stock_thm4;
	end
	close(h)
	%%
	%% 	Now i compute some statistics (mean, median, std, deciles, HPD intervals)
	%%
	MeanMean = zeros(nvar,1);
	MedianMean = zeros(nvar,1);
	StdMean = zeros(nvar,1);
	DistribMean = zeros(nvar,9);
	HPDMean = zeros(nvar,2);
	tmp = zeros(B,1);
	for i = 1:nvar
	  	StartLine = 0;
	  	for file = 1:sfil_thm1 
	  	  	instr = [M_.fname '_thm1' int2str(file)];
	  	  	eval(['load ' instr]);
	  	  	DeProfundis = size(stock_thm1,2);
	  	  	tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm1(i,:));
	  	  	StartLine = StartLine+DeProfundis;
	  	end
	  	tmp = sort(tmp);
	  	MeanMean(i) = mean(tmp);
	  	MedianMean(i) = tmp(round(B*0.5));
	  	StdMean(i) = std(tmp);
	  	DistribMean(i,:) = reshape(tmp(deciles),1,9);
	  	tt = floor(options_.mh_conf_sig*B);
	  	a = 1; 
	  	b = tt;
	  	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	  	while b <= B
	    	tmp1 = [a;b;tmp(b)-tmp(a)];
	    	a = a + 1;
	    	b = b + 1;
			if tmp1(3,1) < tmp2(3,1)
				tmp2 = tmp1;
			end
		end
		HPDMean(i,1) = tmp(tmp2(1,1));
		HPDMean(i,2) = tmp(tmp2(2,1));
	end
	disp(' ')
	disp(' ')
	disp('POSTERIOR THEORETICAL EXPECTATION')
	disp(' ')
	titre = sprintf('%15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
		'Variables',...
		'mean  ',...
		'median',...
		'std   ',...
		'HPDinf',...
		'HPDsup');
	disp(titre)
	for i=1:nvar
		disp(sprintf('%15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
			 deblank(M_.endo_names(ivar(i),:)), ...
			 MeanMean(i),...
			 MedianMean(i),...
			 StdMean(i),...
			 HPDMean(i,1),...
			 HPDMean(i,2)));
		eval(['oo_.PosteriorTheoreticalMoment.Expectation.Mean.' deblank(M_.endo_names(ivar(i),:)) ' = MeanMean(i);']);
		eval(['oo_.PosteriorTheoreticalMoment.Expectation.Median.' deblank(M_.endo_names(ivar(i),:)) ' = MedianMean(i);']);
		eval(['oo_.PosteriorTheoreticalMoment.Expectation.Std.' deblank(M_.endo_names(ivar(i),:)) ' = StdMean(i);']);
		eval(['oo_.PosteriorTheoreticalMoment.Expectation.HPDinf.' deblank(M_.endo_names(ivar(i),:)) ' = HPDMean(i,1);']);
		eval(['oo_.PosteriorTheoreticalMoment.Expectation.HPDsup.' deblank(M_.endo_names(ivar(i),:)) ' = HPDMean(i,2);']);
	end
	if TeX
		fidTeX = fopen([M_.fname '_PosteriorTheoreticalExpectation.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{l|ccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,' Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvar
	    	fprintf(fidTeX,' $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		    	deblank(M_.endo_names_tex(ivar(i),:)), ...
				MeanMean(i),...
				MedianMean(i),...
				StdMean(i),...
				HPDMean(i,1),...
				HPDMean(i,2));
		end
	    fprintf(fidTeX,'\\hline\\hline \n');
	    fprintf(fidTeX,'\\end{tabular}\n ');    
	    fprintf(fidTeX,'\\caption{Posterior theoretical expectation.}\n ');
	    fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalExpectation}\n');
	    fprintf(fidTeX,'\\end{table}\n');
	    fprintf(fidTeX,'} \n');
	    fprintf(fidTeX,'%% End of TeX file.\n');
	    fclose(fidTeX);
	end	
	MeanVariance = zeros(nvar,nvar,1);
	MedianVariance = zeros(nvar,nvar,1);
	StdVariance = zeros(nvar,nvar,1);
	DistribVariance = zeros(nvar,nvar,9);
	HPDVariance = zeros(nvar,nvar,2);
	for i = 1:nvar
		for j=1:nvar
	    	StartLine = 0;
	    	tmp = zeros(B,1);
	    	for file = 1:sfil_thm2 
				instr = [M_.fname '_thm2' int2str(file)];
				eval(['load ' instr]);
				DeProfundis = size(stock_thm2,3);
				tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,j,:));
				StartLine = StartLine+DeProfundis;
			end
			tmp = sort(tmp);
			MeanVariance(i,j) = mean(tmp);
			MedianVariance(i,j) = tmp(round(B*0.5));
			StdVariance(i,j) = std(tmp);
			DistribVariance(i,j,:) = reshape(tmp(deciles),1,1,9);
			tt = floor(options_.mh_conf_sig*B);
			a = 1; 
			b = tt;
			tmp2 = [1;tt;tmp(tt)-tmp(1)];
			while b <= B
				tmp1 = [a;b;tmp(b)-tmp(a)];
				a = a + 1;
				b = b + 1;
				if tmp1(3,1) < tmp2(3,1)
					tmp2 = tmp1;     
				end    
			end
	      	HPDVariance(i,j,1) = tmp(tmp2(1,1));
	      	HPDVariance(i,j,2) = tmp(tmp2(2,1));
		end	
	end
	disp(' ')
	disp(' ')
	disp('POSTERIOR THEORETICAL VARIANCES AND COVARIANCES')
	disp(' ')
	titre = sprintf('%15s \t %15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
		'Variables',...
		'Variables',...
		'mean  ',...
		'median',...
		'std   ',...
		'HPDinf',...
		'HDPsup');
	disp(titre)
	for i=1:nvar
		for j=i:nvar
			disp(sprintf('%15s \t %15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
				deblank(M_.endo_names(ivar(i),:)), ...
				deblank(M_.endo_names(ivar(j),:)), ...
				MeanVariance(i,j),...
				MedianVariance(i,j),...
				StdVariance(i,j),...
				HPDVariance(i,j,1),...
				HPDVariance(i,j,2)));
			eval(['oo_.PosteriorTheoreticalMoment.Variance.Mean.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = MeanVariance(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Variance.Median.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = MedianVariance(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Variance.Std.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = StdVariance(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Variance.HPDinf.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = HPDVariance(i,j,1);']);
			eval(['oo_.PosteriorTheoreticalMoment.Variance.HPDsup.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = HPDVariance(i,j,2);']);
	    end
	end
	if TeX
		fidTeX = fopen([M_.fname '_PosteriorTheoreticalVariance.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,' Variables & Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvar
	    	for j=i:nvar
				fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
					deblank(M_.endo_names_tex(ivar(i),:)), ...
					deblank(M_.endo_names_tex(ivar(j),:)), ...
					MeanVariance(i,j),...
					MedianVariance(i,j),...
					StdVariance(i,j),...
					HPDVariance(i,j,1),...
					HPDVariance(i,j,2));
			end
		end
	    fprintf(fidTeX,'\\hline\\hline \n');
	    fprintf(fidTeX,'\\end{tabular}\n ');    
	    fprintf(fidTeX,'\\caption{Posterior theoretical variances and covariances.}\n ');
	    fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalVariances}\n');
	    fprintf(fidTeX,'\\end{table}\n');
	    fprintf(fidTeX,'} \n');
	    fprintf(fidTeX,'%% End of TeX file.\n');
	    fclose(fidTeX);
	end
	MeanCorrelation = zeros(nvar,nvar,1);
	MedianCorrelation = zeros(nvar,nvar,1);
	StdCorrelation = zeros(nvar,nvar,1);
	DistribCorrelation = zeros(nvar,nvar,9);
	HPDCorrelation = zeros(nvar,nvar,2);
	tmpp	= zeros(B,1);
	tmppp	= zeros(B,1);
	for i = 1:nvar
		for j=1:nvar
	    	StartLine = 0;
	    	tmp = zeros(B,1);
	    	for file = 1:sfil_thm2 
				instr = [M_.fname '_thm2' int2str(file)];
				eval(['load ' instr]);
				DeProfundis = size(stock_thm2,3);
				tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,j,:));
				tmpp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,i,:));
				tmppp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(j,j,:));
				StartLine = StartLine+DeProfundis;
			end
	      	tmp = sort(tmp./sqrt(tmpp.*tmppp));
	      	MeanCorrelation(i,j) = mean(tmp);
	      	MedianCorrelation(i,j) = tmp(round(B*0.5));
	      	StdCorrelation(i,j) = std(tmp);
	      	DistribCorrelation(i,j,:) = reshape(tmp(deciles),1,1,9);
	      	tt = floor(options_.mh_conf_sig*B);
	      	a = 1; 
	      	b = tt;
	      	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	      	while b <= B
				tmp1 = [a;b;tmp(b)-tmp(a)];
				a = a + 1;
				b = b + 1;
				if tmp1(3,1) < tmp2(3,1)
					tmp2 = tmp1;
				end
			end
	      	HPDCorrelation(i,j,1) = tmp(tmp2(1,1));
	      	HPDCorrelation(i,j,2) = tmp(tmp2(2,1));
		end
	end
	clear tmpp tmppp;
	disp(' ')
	disp(' ')
	disp('POSTERIOR THEORETICAL CORRELATIONS')
	disp(' ')
	titre = sprintf('%15s \t %15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
		'Variables',...
		'Variables',...
		'mean  ',...
		'median',...
		'std   ',...
		'HPDinf',...
		'HPDsup');
	disp(titre)
	for i=1:nvar-1
		for j=i+1:nvar
			disp(sprintf('%15s \t %15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
				deblank(M_.endo_names(ivar(i),:)), ...
				deblank(M_.endo_names(ivar(j),:)), ...
	   			MeanCorrelation(i,j),...
	   			MedianCorrelation(i,j),...
	   			StdCorrelation(i,j),...
	   			HPDCorrelation(i,j,1),...
	   			HPDCorrelation(i,j,2)));
			eval(['oo_.PosteriorTheoreticalMoment.Correlation.Mean.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = MeanCorrelation(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Correlation.Median.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = MedianCorrelation(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Correlation.Std.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = StdCorrelation(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Correlation.HPDinf.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = HPDCorrelation(i,j,1);']);
			eval(['oo_.PosteriorTheoreticalMoment.Correlation.HPDsup.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.endo_names(ivar(j),:)) ' = HPDCorrelation(i,j,2);']);
		end
	end
	if TeX
		fidTeX = fopen([M_.fname '_PosteriorTheoreticalCorrelation.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,' Variables & Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvar-1
	    	for j=i+1:nvar
				fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
					deblank(M_.endo_names_tex(ivar(i),:)), ...
					deblank(M_.endo_names_tex(ivar(j),:)), ...
					MeanCorrelation(i,j),...
					MedianCorrelation(i,j),...
					StdCorrelation(i,j),...
					HPDCorrelation(i,j,1),...
					HPDCorrelation(i,j,2));
			end
		end
	    fprintf(fidTeX,'\\hline\\hline \n');
	    fprintf(fidTeX,'\\end{tabular}\n ');    
	    fprintf(fidTeX,'\\caption{Posterior theoretical correlations.}\n ');
	    fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalCorrelations}\n');
	    fprintf(fidTeX,'\\end{table}\n');
	    fprintf(fidTeX,'} \n');
	    fprintf(fidTeX,'%% End of TeX file.\n');
	    fclose(fidTeX);
	end
	MeanDecomp = zeros(nvar,M_.exo_nbr,1);
	MedianDecomp = zeros(nvar,M_.exo_nbr,1);
	StdDecomp = zeros(nvar,M_.exo_nbr,1);
	DistribDecomp = zeros(nvar,M_.exo_nbr,9);
	HPDDecomp = zeros(nvar,M_.exo_nbr,2);
	for i = 1:nvar
		for j=1:M_.exo_nbr
	    	StartLine = 0;
	    	for file = 1:sfil_thm3 
				instr = [M_.fname '_thm3' int2str(file)];
				eval(['load ' instr]);
				DeProfundis = size(stock_thm3,3);
				tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm3(i,j,:));
				StartLine = StartLine+DeProfundis;
			end
	      	tmp = sort(tmp);
	      	MeanDecomp(i,j) = mean(tmp);
	      	MedianDecomp(i,j) = tmp(round(B*0.5));
	      	StdDecomp(i,j) = std(tmp);
	      	DistribDecomp(i,j,:) = reshape(tmp(deciles),1,1,9);
	      	tt = floor(options_.mh_conf_sig*B);
	      	a = 1; 
	      	b = tt;
	      	tmp2 = [1;tt;tmp(tt)-tmp(1)];
			while b <= B
				tmp1 = [a;b;tmp(b)-tmp(a)];
				a = a + 1;
				b = b + 1;
				if tmp1(3,1) < tmp2(3,1)
					tmp2 = tmp1;
				end
			end
			HPDDecomp(i,j,1) = tmp(tmp2(1,1));
			HPDDecomp(i,j,2) = tmp(tmp2(2,1));
		end
	end
	disp(' ')
	disp(' ')
	disp('POSTERIOR THEORETICAL VARIANCE DECOMPOSITION')
	disp(' ')
	titre = sprintf('%15s \t %15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
		'Variables',...
		'Sources',...
		'mean  ',...
		'median',...
		'std   ',...
		'HPDinf',...
		'HDPsup');
	disp(titre)
	for i=1:nvar
		for j=1:M_.exo_nbr
	    	disp(sprintf('%15s \t %15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
				deblank(M_.endo_names(ivar(i),:)), ...
	   			deblank(M_.exo_names(j,:)), ...
	   			MeanDecomp(i,j),...
	   			MedianDecomp(i,j),...
	   			StdDecomp(i,j),...
	   			HPDDecomp(i,j,1),...
	   			HPDDecomp(i,j,2)));
			eval(['oo_.PosteriorTheoreticalMoment.Decomp.Mean.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.exo_names(j,:)) ' = MeanDecomp(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Decomp.Median.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.exo_names(j,:)) ' = MedianDecomp(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Decomp.Std.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.exo_names(j,:)) ' = StdDecomp(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.Decomp.HPDinf.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.exo_names(j,:)) ' = HPDDecomp(i,j,1);']);
			eval(['oo_.PosteriorTheoreticalMoment.Decomp.HPDsup.' deblank(M_.endo_names(ivar(i),:)) '_' deblank(M_.exo_names(j,:)) ' = HPDDecomp(i,j,2);']);
		end
	end
	if TeX
		fidTeX = fopen([M_.fname '_PosteriorTheoreticalVarianceDecomposition.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,' Variables & Sources & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvar
	    	for j=1:M_.exo_nbr
				fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
					deblank(M_.endo_names_tex(ivar(i),:)), ...
					deblank(M_.exo_names_tex(j,:)), ...
					MeanDecomp(i,j),...
					MedianDecomp(i,j),...
					StdDecomp(i,j),...
					HPDDecomp(i,j,1),...
					HPDDecomp(i,j,2));
			end
	    end
	    fprintf(fidTeX,'\\hline\\hline \n');
	    fprintf(fidTeX,'\\end{tabular}\n ');    
	    fprintf(fidTeX,'\\caption{Posterior theoretical variance decomposition.}\n ');
	    fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalVarianceDecomposition}\n');
	    fprintf(fidTeX,'\\end{table}\n');
	    fprintf(fidTeX,'} \n');
	    fprintf(fidTeX,'%% End of TeX file.\n');
	    fclose(fidTeX);
	end
	MeanAutoCorr = zeros(nvar,nar,1);
	MedianAutoCorr = zeros(nvar,nar,1);
	StdAutoCorr = zeros(nvar,nar,1);
	DistribAutoCorr = zeros(nvar,nar,9);
	HPDAutoCorr = zeros(nvar,nar,2);
	for i = 1:nvar
		for j=1:nar
	    	StartLine = 0;
			for file = 1:sfil_thm4 
				instr = [M_.fname '_thm4' int2str(file)];
				eval(['load ' instr]);
				DeProfundis = size(stock_thm4,3);
				tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm4(i,j,:));
				StartLine = StartLine+DeProfundis;
			end
			tmp = sort(tmp);
	      	MeanAutoCorr(i,j) = mean(tmp);
	      	MedianAutoCorr(i,j) = tmp(round(B*0.5));
	      	StdAutoCorr(i,j) = std(tmp);
	      	DistribAutoCorr(i,j,:) = reshape(tmp(deciles),1,1,9);
	      	tt = floor(options_.mh_conf_sig*B);
	      	a = 1; 
	      	b = tt;
	      	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	      	while b <= B
				tmp1 = [a;b;tmp(b)-tmp(a)];
				a = a + 1;
				b = b + 1;
				if tmp1(3,1) < tmp2(3,1)
					tmp2 = tmp1;
				end
			end
	      	HPDAutoCorr(i,j,1) = tmp(tmp2(1,1));
	      	HPDAutoCorr(i,j,2) = tmp(tmp2(2,1));
		end
	end
	disp(' ')
	disp(' ')
	disp('POSTERIOR THEORETICAL AUTOCORRELATION')
	disp(' ')
	titre = sprintf('%15s \t %3s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
		'Variables',...
		'Lag',...
		'mean  ',...
		'median',...
		'std   ',...
		'HPDinf',...
		'HDPsup');
	disp(titre)
	for i=1:nvar
		for j=1:nar
	    	disp(sprintf('%15s \t %3s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
				deblank(M_.endo_names(ivar(i),:)), ...
				[int2str(j) ' '], ...
				MeanAutoCorr(i,j),...
				MedianAutoCorr(i,j),...
				StdAutoCorr(i,j),...
				HPDAutoCorr(i,j,1),...
				HPDAutoCorr(i,j,2)));
			eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Mean.' deblank(M_.endo_names(ivar(i),:)) '_lag' int2str(j) ' = MeanAutoCorr(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Median.' deblank(M_.endo_names(ivar(i),:)) '_lag' int2str(j) ' = MedianAutoCorr(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Std.' deblank(M_.endo_names(ivar(i),:)) '_lag' int2str(j) ' = StdAutoCorr(i,j);']);
			eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.HPDinf.' deblank(M_.endo_names(ivar(i),:)) '_lag' int2str(j) ' = HPDAutoCorr(i,j,1);']);
			eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.HPDsup.' deblank(M_.endo_names(ivar(i),:)) '_lag' int2str(j) ' = HPDAutoCorr(i,j,2);']);
	    end
	end
	if TeX
		fidTeX = fopen([M_.fname '_PosteriorTheoreticalAutocorrelation.TeX'],'w');
		fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
		fprintf(fidTeX,['%% ' datestr(now,0)]);
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,' \n');
		fprintf(fidTeX,'{\\tiny \n');
		fprintf(fidTeX,'\\begin{table}\n');
		fprintf(fidTeX,'\\centering\n');
		fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
		fprintf(fidTeX,'\\hline\\hline \\\\ \n');
		fprintf(fidTeX,' Variables & Lag & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
		fprintf(fidTeX,'\\hline \\\\ \n');
		for i=1:nvar
	    	for j=1:nar
				fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
					deblank(M_.endo_names_tex(ivar(i),:)), ...
					int2str(j), ...
					MeanAutoCorr(i,j),...
					MedianAutoCorr(i,j),...
					StdAutoCorr(i,j),...
					HPDAutoCorr(i,j,1),...
					HPDAutoCorr(i,j,2));
			end
	    end
	    fprintf(fidTeX,'\\hline\\hline \n');
	    fprintf(fidTeX,'\\end{tabular}\n ');    
	    fprintf(fidTeX,'\\caption{Posterior theoretical auto-correlation.}\n ');
	    fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalAutoCorrelation}\n');
	    fprintf(fidTeX,'\\end{table}\n');
	    fprintf(fidTeX,'} \n');
	    fprintf(fidTeX,'%% End of TeX file.\n');
	    fclose(fidTeX);
	end
end% options_.moments_varendo
%% Un dernier petit coup de DsgeLikelihood juste pour remettre les parametres
%% structurels et la matrice de variance-covariance aux valeurs qui
%% correspondent a la moyenne posterieure (en vue d'une utilisation éventuelle
%% de stoch_simul après le Metropolis-Hastings).
[lnprior,cost_flag,ys,trend_coeff] = DsgeLikelihood(post_mean,gend,data);
%% Now I save the seeds (If the user wants to start another MH, he can start from the
%% previous state of the random number generator by using the command "LoadPreviousSeed"
%% before the estimation command)
Seed.NormalDeviates  = randn('state');
Seed.UniformDeviates = rand('state');
save LastSeed Seed;
%% That's All!
%% SA 08-18-2004	* Corrected a bug in forecasts (HPD intervals).
%%					* metropolis.m now displays "true bayesian" smooth shocks. The mean
%%					- across the metropolis draws - of the smooth shocks instead of the 
%%					smooth shocks obtained from the posterior mean are displayed.
%%					* Added "true bayesian" smooth measurement error.
%%					* Added "true bayesian" smooth variables (all the variables in the 
%%					state vector).
%%					* Added deciles for the posterior distribution of the smooth shocks,
%%					variables and measurement errors (green curves).
%% SA 08-19-2004	* Added posterior IRFs.
%% SA 08-21-2004	* Added posterior theoretical moments.
%% SA 08-23-2004   * Added correction to the modified harmonic mean estimator of the
%%                 log-marginal density. The variance of the weighting distribution
%%                 automatically increases if needed (to be revised).)
%% SA 12-02-2004	* Changed architecture for the global structure oo_