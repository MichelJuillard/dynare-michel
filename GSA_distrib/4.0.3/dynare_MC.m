function dynare_MC(var_list_,OutDir)
%
% Adapted by M. Ratto from dynare_estimation.m and posteriorsmoother.m
% (dynare_estimation.m and posteriorsmoother.m are part of DYNARE,
% copyright M. Juillard)
%
% Part of the Sensitivity Analysis Toolbox for DYNARE
%
% Written by Marco Ratto, 2006
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Disclaimer: This software is not subject to copyright protection and is in the public domain. 
% It is an experimental system. The Joint Research Centre of European Commission 
% assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.
%

global M_ options_ oo_ estim_params_ 
global bayestopt_

if options_.filtered_vars ~= 0 & options_.filter_step_ahead == 0
  options_.filter_step_ahead = 1;
end
if options_.filter_step_ahead ~= 0
  options_.nk = max(options_.filter_step_ahead);
else
  options_.nk = 0;
end


nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np ;
npar  = nvx+nvn+ncx+ncn+np;

if isempty(options_.datafile)
  error('ESTIMATION: datafile option is missing')
end

if isempty(options_.varobs)
  error('ESTIMATION: VAROBS is missing')
end

%% Read and demean data 
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);

options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1 & ~options_.logdata
  rawdata = log(rawdata);
end
if options_.prefilter == 1
  bayestopt_.mean_varobs = mean(rawdata,1);
  data = transpose(rawdata-ones(gend,1)*bayestopt_.mean_varobs);
else
  data = transpose(rawdata);
end

if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
	 ' transformation'])
end

offset = npar-np;
fname_=M_.fname;

options_ = set_default_option(options_,'opt_gsa',1);
options_gsa_ = options_.opt_gsa;

if options_gsa_.pprior,
  namfile=[fname_,'_prior'];
else
  namfile=[fname_,'_mc'];
end
load([OutDir,'\',namfile],'lpmat', 'lpmat0', 'istable')
% load(options_.mode_file)
%%
%%
%%
x=[lpmat0(istable,:) lpmat(istable,:)];
clear lpmat lpmat0 istable %iunstable egg yys T
B = size(x,1);
[atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff, aK] = DsgeSmoother(x(1,:)',gend,data);
n1=size(atT,1);

nfil=B/40;
stock_smooth = zeros(M_.endo_nbr,gend,40);
stock_filter = zeros(M_.endo_nbr,gend+1,40);
stock_ys = zeros(40, M_.endo_nbr);
logpo2=zeros(B,1);
%%
h = waitbar(0,'MC smoother ...');
delete([OutDir,'\',namfile,'_*.mat'])
ib=0;
ifil=0;
opt_gsa=options_.opt_gsa;
for b=1:B
  ib=ib+1;
  deep = x(b,:)';
  set_all_parameters(deep);
  dr = resol(oo_.steady_state,0);
  %deep(1:offset) = xparam1(1:offset);
  logpo2(b,1) = DsgeLikelihood(deep,gend,data);
  if opt_gsa.lik_only==0,
  [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff, aK] = DsgeSmoother(deep,gend,data);
  stock_smooth(:,:,ib)=atT(1:M_.endo_nbr,:);
  stock_filter(:,:,ib)=filtered_state_vector(1:M_.endo_nbr,:);
  %stock_filter(:,:,ib)=aK(options_.filter_step_ahead,1:M_.endo_nbr,:);
  stock_ys(ib,:)=ys';
  if ib==40,
    ib=0;
    ifil=ifil+1;
    save([OutDir,'\',namfile,'_',num2str(ifil)],'stock_smooth','stock_filter','stock_ys')
    stock_smooth = zeros(M_.endo_nbr,gend,40);
    stock_filter = zeros(M_.endo_nbr,gend+1,40);
    stock_ys = zeros(40, M_.endo_nbr);
  end
  end  
  waitbar(b/B,h,['MC smoother ...',num2str(b),'/',num2str(B)]);
end
close(h)
if opt_gsa.lik_only==0,
if ib>0,
    ifil=ifil+1;
    stock_smooth = stock_smooth(:,:,1:ib);
    stock_filter = stock_filter(:,:,1:ib);
    stock_ys = stock_ys(1:ib,:);
    save([OutDir,'\',namfile,'_',num2str(ifil)],'stock_smooth','stock_filter','stock_ys')
end
end
stock_gend=gend;
stock_data=data;
save([OutDir,'\',namfile],'x','logpo2','stock_gend','stock_data','-append')
