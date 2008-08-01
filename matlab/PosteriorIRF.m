function PosteriorIRF(type)
% Builds posterior IRFs after the MH algorithm. 
% 
% INPUTS 
%   o type       [char]     string specifying the joint density of the
%                           deep parameters ('prior','posterior'). 
%  
% OUTPUTS 
%   None                    (oo_ and plots).
%
% SPECIAL REQUIREMENTS
%   None

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

global options_ estim_params_ oo_ M_ bayestopt_
% Set the number of periods
if isempty(options_.irf) | ~options_.irf 
    options_.irf = 40;
end
% Set varlist if necessary
varlist = options_.varlist;
if isempty(varlist)
    varlist = options_.varobs;
end
options_.varlist = varlist;
nvar = size(varlist,1);
IndxVariables = [];
for i=1:nvar
    idx = strmatch(deblank(varlist(i,:)),M_.endo_names,'exact');
    if isempty(idx)
        disp(['PosteriorIRF :: ' deblank(varlist(i,:)) 'is not a declared endogenous variable!'])
    else
        IndxVariables = [IndxVariables,idx];
    end
end
% Set various parameters & Check or create directories
nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
clear('nvx','nvn','ncx','ncn','np');
nvobs = size(options_.varobs,1);
gend = options_.nobs;
MaxNumberOfPlotPerFigure = 9;
nn = sqrt(MaxNumberOfPlotPerFigure);
MAX_nirfs_dsge = ceil(options_.MaxNumberOfBytes/(options_.irf*nvar*M_.exo_nbr)/8);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    MAX_nirfs_dsgevar = ceil(options_.MaxNumberOfBytes/(options_.irf*nvobs*M_.exo_nbr)/8);
else
    MAX_nirfs_dsgevar = 0;
end
DirectoryName = CheckPath('Output');
if strcmpi(type,'posterior')
  MhDirectoryName = CheckPath('metropolis');
elseif strcmpi(type,'gsa')
  MhDirectoryName = CheckPath('GSA');
else
  MhDirectoryName = CheckPath('prior');
end
if strcmpi(type,'posterior')
  load([ MhDirectoryName '/'  M_.fname '_mh_history.mat'])
  TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
  NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
elseif strcmpi(type,'gsa')
  load([ MhDirectoryName '/'  M_.fname '_prior.mat'],'lpmat0','lpmat','istable')
  x=[lpmat0(istable,:) lpmat(istable,:)];
  clear lpmat istable
  NumberOfDraws=size(x,1);
  B=NumberOfDraws; options_.B = B;
else% type = 'prior'
  NumberOfDraws = 500;
end
if ~strcmpi(type,'gsa')
  B = min([round(.5*NumberOfDraws),500]); options_.B = B;
end
try delete([MhDirectoryName '/' M_.fname '_irf_dsge*.mat'])
catch disp('No _IRFs (dsge) files to be deleted!')
end
try delete([MhDirectoryName '/' M_.fname '_irf_bvardsge*.mat'])
catch disp('No _IRFs (bvar-dsge) files to be deleted!')
end
irun = 0;
IRUN = 0;
irun2 = 0;
NumberOfIRFfiles_dsge = 1;
NumberOfIRFfiles_dsgevar = 1;
ifil2 = 1;
if strcmpi(type,'posterior')
  h = waitbar(0,'Bayesian (posterior) IRFs...');
elseif strcmpi(type,'gsa')
  h = waitbar(0,'GSA (prior) IRFs...');
else
  h = waitbar(0,'Bayesian (prior) IRFs...');
end
% Create arrays
if B <= MAX_nruns
  stock_param = zeros(B, npar);
else
  stock_param = zeros(MAX_nruns, npar);
end
if B >= MAX_nirfs_dsge
  stock_irf_dsge = zeros(options_.irf,nvar,M_.exo_nbr,MAX_nirfs_dsge);
else
  stock_irf_dsge = zeros(options_.irf,nvar,M_.exo_nbr,B);
end
if MAX_nirfs_dsgevar
    if B >= MAX_nirfs_dsgevar
        stock_irf_bvardsge = zeros(options_.irf,nvobs,M_.exo_nbr,MAX_nirfs_dsgevar);
    else
        stock_irf_bvardsge = zeros(options_.irf,nvobs,M_.exo_nbr,B);
    end
    [mYY,mXY,mYX,mXX,Ydata,Xdata] = ...
        var_sample_moments(options_.first_obs,options_.first_obs+options_.nobs-1,...
                           options_.varlag,-1,options_.datafile,options_.varobs);
    NumberOfLags = options_.varlag;
    NumberOfLagsTimesNvobs = NumberOfLags*nvobs;
    if options_.noconstant
        NumberOfParametersPerEquation = NumberOfLagsTimesNvobs;
    else
        NumberOfParametersPerEquation = NumberOfLagsTimesNvobs+1;
    end
    Companion_matrix = diag(ones(nvobs*(NumberOfLags-1),1),-nvobs);
end
b = 0;
nosaddle = 0;
while b<=B
    b = b + 1;
    irun = irun+1;
    irun2 = irun2+1;
    if ~strcmpi(type,'gsa')
        deep = GetOneDraw(type);
    else
        deep = x(b,:);
    end
    stock_param(irun2,:) = deep;  
    set_parameters(deep);
    [dr,info] = resol(oo_.steady_state,0);
    if info(1)
        nosaddle = nosaddle + 1;
        b = b - 1;
        irun = irun-1;
        irun2 = irun2-1;
        if info(1) == 1
            errordef = 'Static variables are not uniquely defined';
        elseif info(1) == 2
            errordef = 'Dll problem';
        elseif info(1) == 3
            errordef = 'No stable trajectory';
        elseif info(1) == 4
            errordef = 'Indeterminacy';
        elseif info(1) == 5
            errordef = 'Rank condition  is not satisfied';
        end
        disp(['PosteriorIRF :: Dynare is unable to solve the model (' errordef ')'])
        continue
    end
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    SS = transpose(chol(SS));
    for i = 1:M_.exo_nbr
        if SS(i,i) > 1e-13
            y=irf(dr,SS(M_.exo_names_orig_ord,i), options_.irf, options_.drop,options_.replic,options_.order);
            if options_.relative_irf
                y = 100*y/cs(i,i);
            end
            for j = 1:nvar
                if max(y(IndxVariables(j),:)) - min(y(IndxVariables(j),:)) > 1e-12 
                    stock_irf_dsge(:,j,i,irun) = transpose(y(IndxVariables(j),:));
                end
            end
        end
    end
    if MAX_nirfs_dsgevar
        IRUN = IRUN+1;
        %tmp_dsgevar = zeros(options_.irf,nvobs*M_.exo_nbr);
        [fval,cost_flag,info,PHI,SIGMAu,iXX] =  DsgeVarLikelihood(deep',gend);
        dsge_prior_weight = M_.params(strmatch('dsge_prior_weight',M_.param_names));
        DSGE_PRIOR_WEIGHT = floor(gend*(1+dsge_prior_weight));
        SIGMA_inv_upper_chol = chol(inv(SIGMAu*gend*(dsge_prior_weight+1))); 
        explosive_var  = 1;
        while explosive_var
            % draw from the marginal posterior of SIGMA
            SIGMAu_draw = rand_inverse_wishart(nvobs, DSGE_PRIOR_WEIGHT-NumberOfParametersPerEquation, ...
                                               SIGMA_inv_upper_chol);
            % draw from the conditional posterior of PHI
            PHI_draw = rand_matrix_normal(NumberOfParametersPerEquation,nvobs, PHI, ...
                                           chol(SIGMAu_draw)', chol(iXX)');
            Companion_matrix(1:nvobs,:) = transpose(PHI_draw(1:NumberOfLagsTimesNvobs,:));
            % Check for stationarity
            explosive_var = any(abs(eig(Companion_matrix))>1.000000001);
        end
        % Get the mean 
% $$$         if options_.noconstant
            mu = zeros(1,nvobs); 
% $$$         else
% $$$             AA = eye(nvobs);
% $$$             for lag=1:NumberOfLags
% $$$                 AA = AA-PHI_draw((lag-1)*nvobs+1:lag*nvobs,:);
% $$$             end
% $$$             mu = transpose(AA\transpose(PHI_draw(end,:)));
% $$$         end
        % Get rotation
        if dsge_prior_weight > 0
            Atheta(oo_.dr.order_var,M_.exo_names_orig_ord) = oo_.dr.ghu*sqrt(M_.Sigma_e);
            A0 = Atheta(bayestopt_.mfys,:);
            [OMEGAstar,SIGMAtr] = qr2(A0');
        end
        SIGMAu_chol = chol(SIGMAu_draw)';
        SIGMAtrOMEGA = SIGMAu_chol*OMEGAstar';
        PHIpower = eye(NumberOfLagsTimesNvobs);
        irfs = zeros (options_.irf,nvobs*M_.exo_nbr);
        tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
        irfs(1,:) = tmp3(:)';
        for t = 2:options_.irf
            PHIpower = Companion_matrix*PHIpower;
            tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
            irfs(t,:)  = tmp3(:)'+kron(ones(1,M_.exo_nbr),mu);
        end
        tmp_dsgevar = kron(ones(options_.irf,1),mu);
        for j = 1:(nvobs*M_.exo_nbr)
            if max(irfs(:,j)) - min(irfs(:,j)) > 1e-10 
                tmp_dsgevar(:,j) = (irfs(:,j));
            end
        end
        if IRUN < MAX_nirfs_dsgevar
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,nvobs,M_.exo_nbr);
        else
            stock_irf_bvardsge(:,:,:,IRUN) = reshape(tmp_dsgevar,options_.irf,nvobs,M_.exo_nbr);
            instr = [MhDirectoryName '/' M_.fname '_irf_bvardsge' ...
                     int2str(NumberOfIRFfiles_dsgevar) '.mat stock_irf_bvardsge;'];,
            eval(['save ' instr]);
            NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1; 
            IRUN =0;
            stock_irf_dsgevar = zeros(options_.irf,nvobs,M_.exo_nbr,MAX_nirfs_dsgevar);
        end
    end
    if irun == MAX_nirfs_dsge | irun == B | b == B
        if b == B
            stock_irf_dsge = stock_irf_dsge(:,:,:,1:irun);
            if MAX_nirfs_dsgevar & (b == B | IRUN == B)
                stock_irf_bvardsge = stock_irf_bvardsge(:,:,:,1:IRUN);
                instr = [MhDirectoryName '/' M_.fname '_irf_bvardsge' ...
                         int2str(NumberOfIRFfiles_dsgevar) ' stock_irf_bvardsge;'];,
                eval(['save ' instr]);
                NumberOfIRFfiles_dsgevar = NumberOfIRFfiles_dsgevar+1;
                irun = 0;
            end
        end
        save([MhDirectoryName '/' M_.fname '_irf_dsge' int2str(NumberOfIRFfiles_dsge) '.mat'],'stock_irf_dsge');
        NumberOfIRFfiles_dsge = NumberOfIRFfiles_dsge+1;
        irun = 0;
    end
    if irun2 == MAX_nruns | b == B
        if b == B
            stock_param = stock_param(1:irun2,:);
        end
        stock = stock_param;
        save([MhDirectoryName '/' M_.fname '_param_irf' int2str(ifil2) '.mat'],'stock');
        ifil2 = ifil2 + 1;
        irun2 = 0;
    end
    waitbar(b/B,h);
end
if nosaddle
   disp(['PosteriorIRF :: Percentage of discarded posterior draws = ' num2str(nosaddle/(B+nosaddle))]) 
end    
close(h);

ReshapeMatFiles('irf_dsge')
if MAX_nirfs_dsgevar
    ReshapeMatFiles('irf_bvardsge')
end

if strcmpi(type,'gsa')
  return
end

IRF_DSGEs = dir([MhDirectoryName '/' M_.fname '_IRF_DSGEs*.mat']);
NumberOfIRFfiles_dsge = length(IRF_DSGEs);

IRF_BVARDSGEs = dir([MhDirectoryName '/' M_.fname '_IRF_BVARDSGEs*.mat']);
NumberOfIRFfiles_dsgevar = length(IRF_BVARDSGEs);



MeanIRF = zeros(options_.irf,nvar,M_.exo_nbr);
MedianIRF = zeros(options_.irf,nvar,M_.exo_nbr);
VarIRF = zeros(options_.irf,nvar,M_.exo_nbr);
DistribIRF = zeros(options_.irf,9,nvar,M_.exo_nbr);
HPDIRF = zeros(options_.irf,2,nvar,M_.exo_nbr);

if options_.TeX
  varlist_TeX = [];
  for i=1:nvar
    varlist_TeX = strvcat(varlist_TeX,M_.endo_names_tex(IndxVariables(i),:));
  end
end

fprintf('MH: Posterior (dsge) IRFs...\n');
tit(M_.exo_names_orig_ord,:) = M_.exo_names;
kdx = 0;

for file = 1:NumberOfIRFfiles_dsge
  load([MhDirectoryName '/' M_.fname '_IRF_DSGEs' int2str(file) '.mat']);
  for i = 1:M_.exo_nbr
    for j = 1:nvar
        for k = 1:size(STOCK_IRF_DSGE,1)
            kk = k+kdx;
            [MeanIRF(kk,j,i),MedianIRF(kk,j,i),VarIRF(kk,j,i),HPDIRF(kk,:,j,i),...
             DistribIRF(kk,:,j,i)] = posterior_moments(squeeze(STOCK_IRF_DSGE(k,j,i,:)),0,options_.mh_conf_sig);
        end
    end
  end
  kdx = kdx + size(STOCK_IRF_DSGE,1);
end

clear STOCK_IRF_DSGE;

for i = 1:M_.exo_nbr
  for j = 1:nvar
    name = [deblank(M_.endo_names(IndxVariables(j),:)) '_' deblank(tit(i,:))];
    eval(['oo_.PosteriorIRF.dsge.Mean.' name ' = MeanIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Median.' name ' = MedianIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Var.' name ' = VarIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Distribution.' name ' = DistribIRF(:,:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.HPDinf.' name ' = HPDIRF(:,1,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.HPDsup.' name ' = HPDIRF(:,2,j,i);']);
  end
end


if MAX_nirfs_dsgevar
    MeanIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    MedianIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    VarIRFdsgevar = zeros(options_.irf,nvar,M_.exo_nbr);
    DistribIRFdsgevar = zeros(options_.irf,9,nvar,M_.exo_nbr);
    HPDIRFdsgevar = zeros(options_.irf,2,nvar,M_.exo_nbr);    
    fprintf('MH: Posterior (bvar-dsge) IRFs...\n');
    tit(M_.exo_names_orig_ord,:) = M_.exo_names;
    kdx = 0;
    for file = 1:NumberOfIRFfiles_dsgevar
        load([MhDirectoryName '/' M_.fname '_IRF_BVARDSGEs' int2str(file) '.mat']);
        for i = 1:M_.exo_nbr
            for j = 1:nvar
                for k = 1:size(STOCK_IRF_BVARDSGE,1)
                    kk = k+kdx;
                    [MeanIRFdsgevar(kk,j,i),MedianIRFdsgevar(kk,j,i),VarIRFdsgevar(kk,j,i),...
                     HPDIRFdsgevar(kk,:,j,i),DistribIRFdsgevar(kk,:,j,i)] = ...
                        posterior_moments(squeeze(STOCK_IRF_BVARDSGE(k,j,i,:)),0,options_.mh_conf_sig);
                end
            end
        end
        kdx = kdx + size(STOCK_IRF_BVARDSGE,1);
    end
    clear STOCK_IRF_BVARDSGE; 
    for i = 1:M_.exo_nbr
        for j = 1:nvar
            name = [deblank(M_.endo_names(IndxVariables(j),:)) '_' deblank(tit(i,:))];
            eval(['oo_.PosteriorIRF.bvardsge.Mean.' name ' = MeanIRFdsgevar(:,j,i);']);
            eval(['oo_.PosteriorIRF.bvardsge.Median.' name ' = MedianIRFdsgevar(:,j,i);']);
            eval(['oo_.PosteriorIRF.bvardsge.Var.' name ' = VarIRFdsgevar(:,j,i);']);
            eval(['oo_.PosteriorIRF.bvardsge.Distribution.' name ' = DistribIRFdsgevar(:,:,j,i);']);
            eval(['oo_.PosteriorIRF.bvardsge.HPDinf.' name ' = HPDIRFdsgevar(:,1,j,i);']);
            eval(['oo_.PosteriorIRF.bvardsge.HPDsup.' name ' = HPDIRFdsgevar(:,2,j,i);']);
        end
    end
end
%%
%% 	Finally I build the plots.
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
      if ~MAX_nirfs_dsgevar
          h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
          set(h1,'FaceColor',[.9 .9 .9]);
          set(h1,'BaseValue',min(HPDIRF(:,1,j,i)));
          hold on
          h2 = area(1:options_.irf,HPDIRF(:,1,j,i),'FaceColor',[1 1 1],'BaseValue',min(HPDIRF(:,1,j,i)));
          set(h2,'FaceColor',[1 1 1]);
          set(h2,'BaseValue',min(HPDIRF(:,1,j,i)));
          plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
          % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);          
          box on
          axis tight
          xlim([1 options_.irf]);
          hold off
      else    
          h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
          set(h1,'FaceColor',[.9 .9 .9]);
          set(h1,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
          hold on;
          h2 = area(1:options_.irf,HPDIRF(:,1,j,i));
          set(h2,'FaceColor',[1 1 1]);
          set(h2,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
          plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
          % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
          plot(1:options_.irf,MeanIRFdsgevar(:,j,i),'--k','linewidth',2)
          plot(1:options_.irf,HPDIRFdsgevar(:,1,j,i),'--k','linewidth',1)
          plot(1:options_.irf,HPDIRFdsgevar(:,2,j,i),'--k','linewidth',1)
          box on
          axis tight
          xlim([1 options_.irf]);
          hold off
      end
      name = deblank(varlist(j,:));
      NAMES = strvcat(NAMES,name);
      if options_.TeX
        texname = deblank(varlist_TeX(j,:));
        TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
      end
      title(name,'Interpreter','none')
    end
    if subplotnum == MaxNumberOfPlotPerFigure | (j == nvar  & subplotnum> 0)
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