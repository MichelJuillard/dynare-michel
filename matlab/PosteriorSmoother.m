function PosteriorSmoother(Y,gend)
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
ifil1 = 1;
ifil2 = 1;
ifil3 = 1;
ifil4 = 1;
h = waitbar(0,'Bayesian smoother...');
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
    stock_filter = zeros(naK,endo_nbr,gend+1,B);
  else
    stock_filter = zeros(naK,endo_nbr,gend+1,MAX_naK);
  end
end
for b=1:B
  deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
  set_all_parameters(deep);
  dr = resol(oo_.steady_state,0);
  [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(deep,gend,Y)
  
  stock_smooth(:,:,b) = alphahat;
  if nvx
    stock_innov(:,:,b)  = etahat;
  end
  if nvn
    stock_error(:,:,b)  = epsilonhat;
  end
  if naK
    stock_filter(:,:,:,b) = aK;
  end

  irun1 = irun1 + 1;
  irun2 = irun2 + 1;
  irun3 = irun3 + 1;
  irun4 = irun4 + 1;

  if irun1 > MAX_nsmoo | b == B
    if b == B
      stock_smooth = stock_smoo(:,:,1:irun1);
    end
    stock = stock_smooth;
    save([DirectoryName M_.fname '_smooth' int2str(ifil1)],'stock');
    ifil1 = ifil1 + 1;
    irun1 = 1;
  end
  
  if nvx & (irun2 > MAX_inno | b == B)
    if b == B
      stock_innov = stock_innov(:,:,1:irun2);
    end
    stock = stock_inno'
    save([DirectoryName M_.fname '_inno' int2str(ifil2)],'stock');
    ifil2 = ifil2 + 1;
    irun2 = 1;
  end
    
  if nvn & (irun3 > MAX_error | b == B)
    if b == B
      stock_error = stock_error(:,:,1:irun3);
    end
    stock = stock_error;
    save([DirectoryName M_.fname '_error' int2str(ifil3)],'stock');
    ifil3 = ifil3 + 1;
    irun3 = 1;
  end
    
  if naK & (irun3 > MAX_naK | b == B)
    if b == B
      stock_filter = stock_filter(:,:,:,1:irun4);
    end
    stock = stock_filter;
    save([DirectoryName M_.fname '_filter' int2str(ifil4)],'stock');
    ifil4 = ifil4 + 1;
    irun4 = 1;
  end
    
  waitbar(b/B,h);
end
close(h)

