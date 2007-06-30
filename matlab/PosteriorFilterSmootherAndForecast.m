function PosteriorFilterSmootherAndForecast(Y,gend, type)
% stephane.adjemian@ens.fr [09-25-2005]
global options_ estim_params_ oo_ M_ bayestopt_

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
iendo = 1:endo_nbr;
i_last_obs = gend+(1-M_.maximum_endo_lag:0);
horizon = options_.forecast;
maxlag = M_.maximum_endo_lag;
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
B = min(1200, round(0.25*NumberOfDraws));
%%
MAX_nruns = min(B,ceil(options_.MaxNumberOfBytes/(npar+2)/8));
MAX_nsmoo = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8));
MAX_ninno = min(B,ceil(MaxNumberOfBytes/(exo_nbr*gend)/8));
MAX_nerro = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8));
if naK
  MAX_naK   = min(B,ceil(MaxNumberOfBytes/(size(options_.varobs,1)* ...
					   length(options_.filter_step_ahead)*gend)/8));
end
if horizon
  MAX_nforc1 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/8));
  MAX_nforc2 = min(B,ceil(MaxNumberOfBytes/((endo_nbr)*(horizon+maxlag))/ ...
			  8));
  IdObs    = bayestopt_.mfys;

end

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
irun6 = 1;
irun7 = 1;
ifil1 = 1;
ifil2 = 1;
ifil3 = 1;
ifil4 = 1;
ifil5 = 1;
ifil6 = 1;
ifil7 = 1;

h = waitbar(0,'Bayesian smoother...');

stock_param = zeros(MAX_nruns, npar);
stock_logpo = zeros(MAX_nruns,1);
stock_ys = zeros(MAX_nruns,endo_nbr);
if options_.smoother
  stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
  stock_innov  = zeros(exo_nbr,gend,B);
  stock_error = zeros(nvobs,gend,MAX_nerro);
  if options_.filter_step_ahead
    stock_filter = zeros(naK,endo_nbr,gend+options_.filter_step_ahead(end),MAX_naK);
  end
  if options_.forecast
    stock_forcst_mean = zeros(endo_nbr,horizon+maxlag,MAX_forcst1);
    stock_forcst_total = zeros(endo_nbr,horizon+maxlag,MAX_forcst2);
  end
end
for b=1:B
  %deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
  [deep, logpo] = GetOneDraw(type);
  set_all_parameters(deep);
  dr = resol(oo_.steady_state,0);
  [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = ...
      DsgeSmoother(deep,gend,Y);
  
  if options_.loglinear
    stock_smooth(dr.order_var,:,irun1) = alphahat(1:endo_nbr,:)+ ...
	repmat(log(dr.ys(dr.order_var)),1,gend);
  else
    stock_smooth(dr.order_var,:,irun1) = alphahat(1:endo_nbr,:)+ ...
	repmat(dr.ys(dr.order_var),1,gend);
  end    
  if nvx
    stock_innov(:,:,irun2)  = etahat;
  end
  if nvn
    stock_error(:,:,irun3)  = epsilonhat;
  end
  if naK
    stock_filter(:,dr.order_var,:,irun4) = aK(options_.filter_step_ahead,1:endo_nbr,:);
  end
  stock_param(irun5,:) = deep;
  stock_logpo(irun5,1) = logpo;
  stock_ys(irun5,:) = SteadyState';

  if horizon
    yyyy = alphahat(iendo,i_last_obs);
    yf = forcst2a(yyyy,dr,zeros(horizon,exo_nbr));
    if options_.prefilter == 1
      yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
				       horizon+maxlag,1);
    end
    yf(:,IdObs) = yf(:,IdObs)+(gend+[1-maxlag:horizon]')*trend_coeff';
    if options_.loglinear == 1
      yf = yf+repmat(log(SteadyState'),horizon+maxlag,1);
      yf = exp(yf);
    else
      yf = yf+repmat(SteadyState',horizon+maxlag,1);
    end
    yf1 = forcst2(yyyy,horizon,dr,1);
    if options_.prefilter == 1
      yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
	  repmat(bayestopt_.mean_varobs',[horizon+maxlag,1,1]);
    end
    yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-maxlag:horizon]')* ...
					   trend_coeff',[1,1,1]);
    if options_.loglinear == 1
      yf1 = yf1 + repmat(log(SteadyState'),[horizon+maxlag,1,1]);
      yf1 = exp(yf1);
    else
      yf1 = yf1 + repmat(SteadyState',[horizon+maxlag,1,1]);
    end

    stock_forcst_mean(:,:,irun6) = yf;
    stock_forcst_total(:,:,irun7) = yf1;
  end
  
  irun1 = irun1 + 1;
  irun2 = irun2 + 1;
  irun3 = irun3 + 1;
  irun4 = irun4 + 1;
  irun5 = irun5 + 1;
  irun6 = irun6 + 1;
  irun7 = irun7 + 1;

  if irun1 > MAX_nsmoo
    stock = stock_smooth;
    save([DirectoryName '/' M_.fname '_smooth' int2str(ifil1)],'stock');
    ifil1 = ifil1 + 1;
    irun1 = 1;
  end
  
  if nvx & irun2 > MAX_ninno
    stock = stock_innov;
    save([DirectoryName '/' M_.fname '_inno' int2str(ifil2)],'stock');
    ifil2 = ifil2 + 1;
    irun2 = 1;
  end
    
  if nvn & irun3 > MAX_error
    stock = stock_error;
    save([DirectoryName '/' M_.fname '_error' int2str(ifil3)],'stock');
    ifil3 = ifil3 + 1;
    irun3 = 1;
  end
    
  if naK & irun4 > MAX_naK
    stock = stock_filter;
    save([DirectoryName '/' M_.fname '_filter' int2str(ifil4)],'stock');
    ifil4 = ifil4 + 1;
    irun4 = 1;
  end
    
  if irun5 > MAX_nruns
    stock = stock_param;
    save([DirectoryName '/' M_.fname '_param' int2str(ifil5)],'stock','stock_logpo','stock_ys');
    ifil5 = ifil5 + 1;
    irun5 = 1;
  end

  if irun6 > MAX_nforc1
    stock = stock_forcst_mean;
    save([DirectoryName '/' M_.fname '_forc_mean' int2str(ifil6)],'stock');
    ifil6 = ifil6 + 1;
    irun6 = 1;
  end

  if irun7 > MAX_nforc2
    stock = stock_forcst_total;
    save([DirectoryName '/' M_.fname '_forc_total' int2str(ifil7)],'stock');
    ifil6 = ifil6 + 1;
    irun6 = 1;
  end

  waitbar(b/B,h);
end
close(h)

stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data'],'stock_gend','stock_data');

ifil1=25;
pm3(endo_nbr,gend,ifil1,B,'Smoothed variables',...
    M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
    'names2','smooth',[M_.fname '/metropolis'],'_smooth')