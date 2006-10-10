function posterior_irf(type)
% Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o type       [char]     string specifying the joint density of the
%   deep parameters ('prior','posterior'). 
%  
% OUTPUTS 
%   None (oo_ and plots).
%
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
global options_ estim_params_ oo_ M_ dsge_prior_weight
nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
%
MaxNumberOfPlotPerFigure = 9;% The square root must be an integer!
nn = sqrt(MaxNumberOfPlotPerFigure);
DirectoryName = CheckPath('Output');
if strcmpi(type,'posterior')
  MhDirectoryName = CheckPath('metropolis');
else
  MhDirectoryName = CheckPath('prior');
end  
MAX_nirfs = ceil(options_.MaxNumberOfBytes/(options_.irf*length(oo_.steady_state)*M_.exo_nbr)/8)+50;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

if strcmpi(type,'posterior')
  load([ MhDirectoryName '/'  M_.fname '_mh_history'])
  TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
  NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
else% type = 'prior'
  NumberOfDraws = 500;
end
B = min([round(.5*NumberOfDraws),500]); options_.B = B;
try delete([MhDirectoryName '\' M_.fname '_IRFs*']);
catch disp('No _IRFs files to be deleted!')
end
irun = 0;
irun2 = 0;
NumberOfIRFfiles = 1;
ifil2 = 1;
if strcmpi(type,'posterior')
  h = waitbar(0,'Bayesian (posterior) IRFs...');
else
  h = waitbar(0,'Bayesian (prior) IRFs...');
end
if B <= MAX_nruns
  stock_param = zeros(B, npar);
else
  stock_param = zeros(MAX_nruns, npar);
end
if B >= MAX_nirfs
  stock_irf = zeros(options_.irf,M_.endo_nbr,M_.exo_nbr,MAX_nirfs);
else
  stock_irf = zeros(options_.irf,M_.endo_nbr,M_.exo_nbr,B);
end
for b=1:B
  irun = irun+1;
  irun2 = irun2+1;
  deep = GetOneDraw(type);
  stock_param(irun2,:) = deep;  
  set_parameters(deep);
  dr = resol(oo_.steady_state,0);
  SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
  SS = transpose(chol(SS));
  for i = 1:M_.exo_nbr
    if SS(i,i) > 1e-13
      y=irf(dr,SS(M_.exo_names_orig_ord,i), options_.irf, options_.drop,options_.replic,options_.order);
      if options_.relative_irf
        y = 100*y/cs(i,i);
      end
      for j = 1:M_.endo_nbr
        if max(y(j,:)) - min(y(j,:)) > 1e-10 
          stock_irf(:,j,i,irun) = transpose(y(j,:));
        end
      end
    end
  end
  if irun == MAX_nirfs | irun == B | b == B
    if b == B
      stock_irf = stock_irf(:,:,:,1:irun);
    end
    save([MhDirectoryName '/' M_.fname '_irf' int2str(NumberOfIRFfiles)],'stock_irf');
    NumberOfIRFfiles = NumberOfIRFfiles+1;
    irun = 0;
  end
  if irun2 > MAX_nruns | b == B
    if b == B
      stock_param = stock_param(1:irun2,:);
    end
    stock = stock_param;
    save([MhDirectoryName '/' M_.fname '_param_irf' int2str(ifil2)],'stock');
    ifil2 = ifil2 + 1;
    irun2 = 1;
  end
  waitbar(b/B,h);
end
NumberOfIRFfiles = NumberOfIRFfiles-1;
ifil2 = ifil2-1;
close(h);

ReshapeMatFiles('irf')

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
kdx = 0;
for file = 1:NumberOfIRFfiles
  load([MhDirectoryName '/' M_.fname '_IRFs' int2str(file)]);
  for i = 1:M_.exo_nbr
    for j = 1:nvar
      for k = 1:size(STOCK_IRF,1)
        kk = k+kdx;
        [MeanIRF(kk,j,i),MedianIRF(kk,j,i),VarIRF(kk,j,i),HPDIRF(kk,:,j,i),DistribIRF(kk,:,j,i)] = ...
          posterior_moments(squeeze(STOCK_IRF(k,SelecVariables(j),i,:)),0);
      end
    end
  end
  kdx = kdx + size(STOCK_IRF,1);
end
clear STOCK_IRF;

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
  fidTeX = fopen([DirectoryName '/' M_.fname '_BayesianIRF.TeX'],'w');
  fprintf(fidTeX,'%% TeX eps-loader file generated by PosteriorIRF.m (Dynare).\n');
  fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
  fprintf(fidTeX,' \n');
  titTeX(M_.exo_names_orig_ord,:) = M_.exo_names_tex;
end
%%
subplotnum = 0;
for i=1:M_.exo_nbr
  NAMES = [];
  if options_.TeX
    TEXNAMES = [];
  end
  figunumber = 0;
  for j=1:nvar
    if max(abs(MeanIRF(:,j,i))) > 10^(-6)
      subplotnum = subplotnum+1;
      if options_.nograph
        if subplotnum == 1 & options_.relative_irf
          hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)],'Visible','off');
        elseif subplotnum == 1 & ~options_.relative_irf
          hh = figure('Name',['Orthogonalized shock to ' tit(i,:)],'Visible','off');
        end
      else
        if subplotnum == 1 & options_.relative_irf
          hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)]);
        elseif subplotnum == 1 & ~options_.relative_irf
          hh = figure('Name',['Orthogonalized shock to ' tit(i,:)]);
        end
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
    if subplotnum == MaxNumberOfPlotPerFigure | (j == nvar  & subplotnum>0)
      figunumber = figunumber+1;
      set(hh,'visible','on')
      eval(['print -depsc2 ' DirectoryName '/'  M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber)]);
      eval(['print -dpdf ' DirectoryName '/' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber)]);
      saveas(hh,[DirectoryName '/' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:))  '_' int2str(figunumber) '.fig']);
      set(hh,'visible','off')
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
    end
  end% loop over selected endo_var
end% loop over exo_var
%%
if options_.TeX
  fprintf(fidTeX,'%% End of TeX file.\n');
  fclose(fidTeX);
end
fprintf('MH: Posterior IRFs, done!\n');
