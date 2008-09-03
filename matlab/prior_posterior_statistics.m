function prior_posterior_statistics(type,Y,gend)

% function PosteriorFilterSmootherAndForecast(Y,gend, type)
% Computes posterior filter smoother and forecasts
%
% INPUTS
%    type:   posterior
%            prior
%            gsa
%    Y:      data
%    gend:   number of observations 
%    
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2008 Dynare Team
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

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
naK = length(options_.filter_step_ahead);
%%
MaxNumberOfBytes=options_.MaxNumberOfBytes;
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
nvobs     = size(options_.varobs,1);
iendo = 1:endo_nbr;
horizon = options_.forecast;
% moments_varendo = options_.moments_varendo;
filtered_vars = options_.filtered_vars;
if horizon
    i_last_obs = gend+(1-M_.maximum_endo_lag:0);
end
maxlag = M_.maximum_endo_lag;
%%
DirectoryName = CheckPath('metropolis');
load([ DirectoryName '/'  M_.fname '_mh_history'])
FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; 
TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles; 
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
clear record;
if ~isempty(options_.subdraws)
    B = options_.subdraws;
    if B > NumberOfDraws
        B = NumberOfDraws;
    end
else
    B = min(1200, round(0.25*NumberOfDraws));
end

%%
MAX_nruns = min(B,ceil(MaxNumberOfBytes/(npar+2)/8));
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
MAX_momentsno = min(B,ceil(MaxNumberOfBytes/(get_moments_size(options_)*8)));
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

irun = ones(7,1);
ifil = zeros(7,1);

if exist('OCTAVE_VERSION')
    diary off;
else
    h = waitbar(0,'Taking subdraws...');
end

stock_param = zeros(MAX_nruns, npar);
stock_logpo = zeros(MAX_nruns,1);
stock_ys = zeros(MAX_nruns,endo_nbr);
run_smoother = 0;
if options_.smoother
  stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
  stock_innov  = zeros(exo_nbr,gend,B);
  stock_error = zeros(nvobs,gend,MAX_nerro);
  stock_filter = zeros(endo_nbr,gend,MAX_nsmoo);
  run_smoother = 1;
end

if options_.filter_step_ahead
    stock_filter_step_ahead = zeros(naK,endo_nbr,gend+ ...
                         options_.filter_step_ahead(end),MAX_naK);
    run_smoother = 1;
end
if options_.forecast
    stock_forcst_mean = zeros(endo_nbr,horizon+maxlag,MAX_nforc1);
    stock_forcst_point = zeros(endo_nbr,horizon+maxlag,MAX_nforc2);
    run_smoother = 1;
end
%if moments_varendo
%    stock_moments = cell(MAX_momentsno,1);
%end
for b=1:B
  [deep, logpo] = GetOneDraw(type);
  set_all_parameters(deep);
  dr = resol(oo_.steady_state,0);
  %if moments_varendo
  %    stock_moments{irun(8)} = compute_model_moments(dr,M_,options_);
  %end
  if run_smoother
      [alphahat,etahat,epsilonhat,alphatilde,SteadyState,trend_coeff,aK] = ...
          DsgeSmoother(deep,gend,Y);
      if options_.loglinear
          stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
              repmat(log(dr.ys(dr.order_var)),1,gend);
          stock_filter(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
              repmat(log(dr.ys(dr.order_var)),1,gend);
      else
          stock_smooth(dr.order_var,:,irun(1)) = alphahat(1:endo_nbr,:)+ ...
              repmat(dr.ys(dr.order_var),1,gend);
          stock_filter(dr.order_var,:,irun(1)) = alphatilde(1:endo_nbr,:)+ ...
              repmat(dr.ys(dr.order_var),1,gend);
      end    
      stock_innov(:,:,irun(2))  = etahat;
      if nvn
          stock_error(:,:,irun(3))  = epsilonhat;
      end
      if naK
          stock_filter_step_ahead(:,dr.order_var,:,irun(4)) = aK(options_.filter_step_ahead,1:endo_nbr,:);
      end

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
          else
              yf1 = yf1 + repmat(SteadyState',[horizon+maxlag,1,1]);
          end

          stock_forcst_mean(:,:,irun(6)) = yf';
          stock_forcst_point(:,:,irun(7)) = yf1';
      end
      
  end
  stock_param(irun(5),:) = deep;
  stock_logpo(irun(5),1) = logpo;
  stock_ys(irun(5),:) = SteadyState';

  irun = irun +  ones(7,1);

  if irun(1) > MAX_nsmoo || b == B
      stock = stock_smooth(:,:,1:irun(1)-1);
      ifil(1) = ifil(1) + 1;
      save([DirectoryName '/' M_.fname '_smooth' int2str(ifil(1)) '.mat'],'stock');
      stock = stock_filter(:,:,1:irun(1)-1);
      save([DirectoryName '/' M_.fname '_filter' int2str(ifil(1)) '.mat'],'stock');
      irun(1) = 1;
  end
  
  if irun(2) > MAX_ninno || b == B
      stock = stock_innov(:,:,1:irun(2)-1);
      ifil(2) = ifil(2) + 1;
      save([DirectoryName '/' M_.fname '_inno' int2str(ifil(2)) '.mat'],'stock');
      irun(2) = 1;
  end
  
  if nvn && (irun(3) > MAX_nerro || b == B)
      stock = stock_error(:,:,1:irun(3)-1);
      ifil(3) = ifil(3) + 1;
      save([DirectoryName '/' M_.fname '_error' int2str(ifil(3)) '.mat'],'stock');
      irun(3) = 1;
  end
  
  if naK && (irun(4) > MAX_naK || b == B)
      stock = stock_filter_step_ahead(:,:,:,1:irun(4)-1);
      ifil(4) = ifil(4) + 1;
      save([DirectoryName '/' M_.fname '_filter_step_ahead' int2str(ifil(4)) '.mat'],'stock');
      irun(4) = 1;
  end
  
  if irun(5) > MAX_nruns || b == B
      stock = stock_param(1:irun(5)-1,:);
      ifil(5) = ifil(5) + 1;
      save([DirectoryName '/' M_.fname '_param' int2str(ifil(5)) '.mat'],'stock','stock_logpo','stock_ys');
      irun(5) = 1;
  end

  if horizon && (irun(6) > MAX_nforc1 || b == B)
      stock = stock_forcst_mean(:,:,1:irun(6)-1);
      ifil(6) = ifil(6) + 1;
      save([DirectoryName '/' M_.fname '_forc_mean' int2str(ifil(6)) '.mat'],'stock');
      irun(6) = 1;
  end

  if horizon && (irun(7) > MAX_nforc2 ||  b == B)
      stock = stock_forcst_point(:,:,1:irun(7)-1);
      ifil(7) = ifil(7) + 1;
      save([DirectoryName '/' M_.fname '_forc_point' int2str(ifil(7)) '.mat'],'stock');
      irun(7) = 1;
  end

  % if moments_varendo && (irun(8) > MAX_momentsno || b == B)
  %    stock = stock_moments(1:irun(8)-1);
  %    ifil(8) = ifil(8) + 1;
  %    save([DirectoryName '/' M_.fname '_moments' int2str(ifil(8)) '.mat'],'stock');
  %    irun(8) = 1;
  % end
  
  if exist('OCTAVE_VERSION')
      printf('Taking subdraws: %3.f%% done\r', b/B);
  else
      waitbar(b/B,h);
  end
end

if exist('OCTAVE_VERSION')
    printf('\n');
    diary on;
else
    close(h)
end

stock_gend=gend;
stock_data=Y;
save([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

if options_.smoother
    pm3(endo_nbr,gend,ifil(1),B,'Smoothed variables',...
	'',M_.endo_names,'tit_tex',M_.endo_names,...
	varlist,'smoothed_variables',[M_.fname '/metropolis'],'_smooth');
    pm3(exo_nbr,gend,ifil(2),B,'Smoothed shocks',...
	'',M_.exo_names,'tit_tex',M_.exo_names,...
	M_.exo_names,'smoothed_shocks',[M_.fname '/metropolis'],'_inno');
    if nvn
        % needs to  be fixed
%        pm3(endo_nbr,gend,ifil(3),B,'Smoothed measurement errors',...
%            M_.endo_names(SelecVariables),M_.endo_names,'tit_tex',M_.endo_names,...
%            'names2','smooth_errors',[M_.fname '/metropolis'],'_error')
    end
end

if options_.filtered_vars
    pm3(endo_nbr,gend,ifil(1),B,'Filtered variables',...
	'',M_.endo_names,'tit_tex',M_.endo_names,...
	varlist,'filtered_current_variables',[M_.fname '/metropolis'], ...
        '_filter');
end

if options_.filter_step_ahead
    pm3(endo_nbr,gend+1,ifil(1),B,'One step ahead forecast',...
	'',M_.endo_names,'tit_tex',M_.endo_names,...
	varlist,'one_step_ahead',[M_.fname '/metropolis'],'_filter_step_ahead');
end

if options_.forecast
    pm3(endo_nbr,horizon+maxlag,ifil(6),B,'Forecasted variables (mean)',...
	'',M_.endo_names,'tit_tex',M_.endo_names,...
	varlist,'mean_forecast',[M_.fname '/metropolis'],'_forc_mean');
    pm3(endo_nbr,horizon+maxlag,ifil(6),B,'Forecasted variables (point)',...
	'',M_.endo_names,'tit_tex',M_.endo_names,...
	varlist,'point_forecast',[M_.fname '/metropolis'],'_forc_point');
end