function PosteriorIRF()
% stephane.adjemian@ens.fr [09-25-2005]
global options_ estim_params_ oo_ M_

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
%%
MaxNumberOfPlotPerFigure = 4;% The square root must be an integer!
nn = sqrt(MaxNumberOfPlotPerFigure);
%%
CheckPath('Plots/IRFs');
CheckPath('metropolis/IRFs');
DirectoryName = CheckPath('metropolis');
load([ DirectoryName '/'  M_.fname '_mh_history'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; 
TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles; 
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
MAX_nirfs = ceil(options_.MaxNumberOfBytes/(options_.irf*length(oo_.steady_state)*M_.exo_nbr)/8)+50;
%%
B = round(0.25*NumberOfDraws);
if B <= MAX_nirfs
  stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,B);
elseif nvn & B > MAX_nirfs
  stock_irf = zeros(options_.irf,size(M_.endo_names,1),M_.exo_nbr,MAX_nirfs);
end
%%
irun = 0;
ifil = 1;
h = waitbar(0,'Bayesian IRFs...');
if B >= MAX_nirfs 
  stock_irf = zeros(options_.irf,M_.endo_nbr,M_.exo_nbr,MAX_nirfs);
else
  stock_irf = zeros(options_.irf,M_.endo_nbr,M_.exo_nbr,B);
end
for b=1:B
  irun = irun+1;
  deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
  M_.params(estim_params_.param_vals(:,1)) = deep(offset+1:end);
  dr = resol(oo_.steady_state,0);
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
  SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
  SS = transpose(chol(SS));
  for i = 1:M_.exo_nbr
    if SS(i,i) > 1e-13
      y=irf(dr,SS(M_.exo_names_orig_ord,i), options_.irf, options_.drop,options_.replic,options_.order);
      if options_.relative_irf
	y = 100*y/cs(i,i);
      end
      for j = 1:M_.endo_nbr%size(M_.endo_names,1)
	if max(y(j,:)) - min(y(j,:)) > 1e-10 
	  stock_irf(:,j,i,irun) = transpose(y(j,:));
	end
      end
    end
  end
  if irun == MAX_nirfs | irun == B | i == B
    if i == B
      stock_irf = stock_irf(:,:,:,1:irun);
    end
    save([DirectoryName '\IRFs\' M_.fname '_irf' int2str(ifil)],'stock_irf');
    ifil = ifil+1;
    irun = 0;
  end
  waitbar(b/B,h);
end
ifil = ifil-1;
close(h)
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
MeanIRF = zeros(options_.irf,nvar,M_.exo_nbr);
MedianIRF = zeros(options_.irf,nvar,M_.exo_nbr);
StdIRF = zeros(options_.irf,nvar,M_.exo_nbr);
DistribIRF = zeros(options_.irf,9,nvar,M_.exo_nbr);
HPDIRF = zeros(options_.irf,2,nvar,M_.exo_nbr);
if options_.TeX
  varlist_TeX = [];
  for i=1:nvar
    varlist_TeX = strvcat(varlist_TeX,M_.endo_names_tex(SelecVariables(i),:));
  end
end
fprintf('MH: Posterior IRFs...\n');
tit(M_.exo_names_orig_ord,:) = M_.exo_names;
for i = 1:M_.exo_nbr
  for j = 1:nvar
    for k = 1:options_.irf
      StartLine = 0;
      tmp = zeros(B,1);
      for file = 1:ifil
	load([DirectoryName '\IRFs\' M_.fname '_irf' int2str(file)]);
	DeProfundis = size(stock_irf,4); 
	tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_irf(k,SelecVariables(j),i,:)); 
	StartLine = StartLine+DeProfundis;
      end
      [MeanIRF(k,j,i),MedianIRF(k,j,i),VarIRF(k,j,i),HPDIRF(k,:,j,i),DistribIRF(k,:,j,i)] = posterior_moments(tmp,0);
    end
    disp(['    Variable: ' deblank(M_.endo_names(SelecVariables(j),:)) ', orthogonalized shock to ' deblank(tit(i,:))])
  end
end  
clear stock_irf;
for i = 1:M_.exo_nbr
  for j = 1:nvar
    name = [deblank(M_.endo_names(SelecVariables(j),:)) '_' deblank(tit(i,:))];
    eval(['oo_.PosteriorIRF.Mean.' name ' = MeanIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.Median.' name ' = MedianIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.Var.' name ' = VarIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.Distribution.' name ' = DistribIRF(:,:,j,i);']);
    eval(['oo_.PosteriorIRF.HPDinf.' name ' = HPDIRF(:,1,j,i);']);
    eval(['oo_.PosteriorIRF.HPDsup.' name ' = HPDIRF(:,2,j,i);']);
  end
end
%%
%% 	Finally i build the plots.
%%
if options_.TeX
  fidTeX = fopen([M_.dname '\Plots\IRFs\' M_.fname '_BayesianIRF.TeX'],'w');
  fprintf(fidTeX,'%% TeX eps-loader file generated by PosteriorIRF.m (Dynare).\n');
  fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
  fprintf(fidTeX,' \n');
  titTeX(M_.exo_names_orig_ord,:) = M_.exo_names_tex;
end
%%
figunumber = 0;
subplotnum = 0;
for i=1:M_.exo_nbr
  NAMES = [];
  if options_.TeX 
    TEXNAMES = []; 
  end
  for j=1:nvar
    if max(abs(MeanIRF(:,j,i))) > 10^(-6)
      subplotnum = subplotnum+1;
      if subplotnum == 1 & options_.relative_irf
	hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)]);
      elseif subplotnum == 1 & ~options_.relative_irf
	hh = figure('Name',['Orthogonalized shock to ' tit(i,:)]);
      end
      set(0,'CurrentFigure',hh)
      subplot(nn,nn,subplotnum);
      plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
      hold on
      for k = 1:9
 	plot(1:options_.irf,DistribIRF(:,k,j,i),'-g','linewidth',0.5)
      end
      plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',1)
      xlim([1 options_.irf]);
      hold off
      name = deblank(varlist(j,:));
      NAMES = strvcat(NAMES,name);
      if options_.TeX
 	texname = deblank(varlist_TeX(j,:));
 	TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
      end
      title(name,'Interpreter','none')
    end
    if subplotnum == MaxNumberOfPlotPerFigure | j == nvar  
      eval(['print -depsc2 ' M_.dname '\Plots\IRFs\'  M_.fname '_Bayesian_IRF_' deblank(tit(i,:))]);
      eval(['print -dpdf ' M_.dname '\Plots\IRFs\' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:))]);
      saveas(hh,[M_.dname '\Plots\IRFs\' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
	fprintf(fidTeX,'\\begin{figure}[H]\n');
	for jj = 1:size(TEXNAMES,1)
	  fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	end    
	fprintf(fidTeX,'\\centering \n');
	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRF_%s}\n',M_.fname,deblank(tit(i,:)));
	if options_.relative_irf
	  fprintf(fidTeX,['\\caption{Bayesian relative IRF.}']);
	else
	  fprintf(fidTeX,'\\caption{Bayesian IRF.}');
	end
	fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s}\n',deblank(tit(i,:)));
	fprintf(fidTeX,'\\end{figure}\n');
	fprintf(fidTeX,' \n');
      end
      subplotnum = 0;
      figunumber = figunumber+1;
    end
  end% loop over selected endo_var
end% loop over exo_var  
%%
if options_.TeX
  fprintf(fidTeX,'%% End of TeX file.\n');
  fclose(fidTeX);
end
fprintf('MH: Posterior IRFs, done!\n');