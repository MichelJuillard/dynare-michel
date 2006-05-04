function PosteriorSmoother(Y,gend, type)
% stephane.adjemian@ens.fr [09-25-2005]
global options_ estim_params_ oo_ M_

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
naK = length(options_.filter_step_ahead);
%%
MaxNumberOfPlotPerFigure = 4;% The square root must be an integer!
MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
nvobs     = size(options_.varobs,1);
nn = sqrt(MaxNumberOfPlotPerFigure);
%%
CheckPath('Plots/');
DirectoryName = CheckPath('metropolis');
load([ DirectoryName '/'  M_.fname '_mh_history'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; 
TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles; 
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
MAX_nsmoo = ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8);
MAX_ninno = ceil(MaxNumberOfBytes/(exo_nbr*gend)/8);
MAX_nerro = ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8);
MAX_naK   = ceil(MaxNumberOfBytes/(size(options_.varobs,1)*length(options_.filter_step_ahead)*gend)/8);
%%
B = round(0.25*NumberOfDraws);
%%
varlist = options_.varlist;
if isempty(varlist)
  varlist = M_.endo_names;
  SelecVariables = transpose(1:M_.endo_nbr);
  nvar = M_.endo_nbr;
else
  nvar = size(varlist,1);
  SelecVariables = [];
  for i=1:nvar
    if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
      SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
    end
  end
end

irun1 = 1;
irun2 = 1;
irun3 = 1;
irun4 = 1;
irun5 = 1;
ifil1 = 1;
ifil2 = 1;
ifil3 = 1;
ifil4 = 1;
ifil5 = 1;
h = waitbar(0,'Bayesian smoother...');
if B <= MAX_nruns
  stock_param = zeros(B, npar);
  stock_logpo = zeros(B,1);
  stock_ys = zeros(B,endo_nbr);
else
  stock_param = zeros(MAX_nruns, npar);
  stock_logpo = zeros(MAX_nruns,1);
  stock_ys = zeros(MAX_nruns,endo_nbr);
end
if options_.smoother
  if B <= MAX_nsmoo
    stock_smooth = zeros(endo_nbr,gend,B);
  else
    stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
  end
  if B <= MAX_ninno 
    stock_innov  = zeros(exo_nbr,gend,B);
  else
    stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
  end
  if nvn & B <= MAX_nerro
    stock_error = zeros(nvobs,gend,B);
  else nvn & B > MAX_nerro
    stock_error = zeros(nvobs,gend,MAX_nerro);
  end
end
if options_.filter_step_ahead ~= 0
  if B <= MAX_naK
    stock_filter = zeros(naK,endo_nbr,gend+options_.filter_step_ahead(end),B);
  else
    stock_filter = zeros(naK,endo_nbr,gend+options_.filter_step_ahead(end),MAX_naK);
  end
end
for b=1:B
  %deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
  [deep, logpo] = GetOneDraw(type);
  set_all_parameters(deep);
  dr = resol(oo_.steady_state,0);
  [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(deep,gend,Y);
  
  stock_smooth(:,:,irun1) = alphahat(1:endo_nbr,:);
  if nvx
    stock_innov(:,:,irun2)  = etahat;
  end
  if nvn
    stock_error(:,:,irun3)  = epsilonhat;
  end
  if naK
    stock_filter(:,:,:,irun4) = aK(options_.filter_step_ahead,1:endo_nbr,:);
  end
  stock_param(irun5,:) = deep;
  stock_logpo(irun5,1) = logpo;
  stock_ys(irun5,:) = SteadyState';

  irun1 = irun1 + 1;
  irun2 = irun2 + 1;
  irun3 = irun3 + 1;
  irun4 = irun4 + 1;
  irun5 = irun5 + 1;

  if irun1 > MAX_nsmoo | b == B
    if b == B
      stock_smooth = stock_smooth(:,:,1:irun1-1);
    end
    stock = stock_smooth;
    save([DirectoryName '/' M_.fname '_smooth' int2str(ifil1)],'stock');
    ifil1 = ifil1 + 1;
    irun1 = 1;
  end
  
  if nvx & (irun2 > MAX_ninno | b == B)
    if b == B
      stock_innov = stock_innov(:,:,1:irun2-1);
    end
    stock = stock_innov;
    save([DirectoryName '/' M_.fname '_inno' int2str(ifil2)],'stock');
    ifil2 = ifil2 + 1;
    irun2 = 1;
  end
    
  if nvn & (irun3 > MAX_error | b == B)
    if b == B
      stock_error = stock_error(:,:,1:irun3-1);
    end
    stock = stock_error;
    save([DirectoryName '/' M_.fname '_error' int2str(ifil3)],'stock');
    ifil3 = ifil3 + 1;
    irun3 = 1;
  end
    
  if naK & (irun4 > MAX_naK | b == B)
    if b == B
      stock_filter = stock_filter(:,:,:,1:irun4-1);
    end
    stock = stock_filter;
    save([DirectoryName '/' M_.fname '_filter' int2str(ifil4)],'stock');
    ifil4 = ifil4 + 1;
    irun4 = 1;
  end
    
  if irun5 > MAX_nruns | b == B
    if b == B
      stock_param = stock_param(1:irun5-1,:);
      stock_logpo = stock_logpo(1:irun5-1,1);
      stock_ys = stock_ys(1:irun5-1,:);
    end
    stock = stock_param;
    save([DirectoryName '/' M_.fname '_param' int2str(ifil5)],'stock','stock_logpo','stock_ys');
    ifil5 = ifil5 + 1;
    irun5 = 1;
  end

  waitbar(b/B,h);
end
close(h)

stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data'],'stock_gend','stock_data');
