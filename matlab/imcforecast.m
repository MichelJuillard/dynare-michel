function icforecast(ptype,cV,cS,cL,H,mcValue,B,ci)
% stephane.adjemian@ens.fr
global options_ oo_ M_
   
xparam = get_posterior_parameters(ptype);
gend   = options_.nobs;

% Read and demean data 
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);
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

set_parameters(xparam);

[atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(xparam,gend,data);

InitState(:,1) = atT(:,end);
[T,R,ys,info] = dynare_resolve;

sQ = sqrt(M_.Sigma_e);

NumberOfStates = length(InitState);
FORCS1 = zeros(NumberOfStates,H+1,B);

for b=1:B
    FORCS1(:,1,b) = InitState;
end

EndoSize = size(M_.endo_names,1);
ExoSize = size(M_.exo_names,1);

n1 = size(cV,1);
n2 = size(cS,1);

if n1 ~= n2
    disp('imcforecast :: Error!')
    disp(['imcforecast :: The number of variables doesn''t match the number of shocks'])
    return
end

idx = [];
jdx = [];

for i = 1:n1
    idx = [idx ; oo_.dr.inv_order_var(strmatch(deblank(cV(i,:)),M_.endo_names,'exact'))];
    jdx = [jdx ; strmatch(deblank(cS(i,:)),M_.exo_names,'exact')];
end
mv = zeros(n1,NumberOfStates);
mu = zeros(ExoSize,n2);
for i=1:n1
    mv(i,idx(i)) = 1;
    mu(jdx(i),i) = 1;
end


if (size(mcValue,2) == 1);
        mcValue = mcValue*ones(1,cL);
else
    cL = size(mcValue,2);
end

randn('state',0);

for b=1:B
    shocks = sQ*randn(ExoSize,H);
    shocks(jdx,:) = zeros(length(jdx),H);
    FORCS1(:,:,b) = mcforecast3(cL,H,mcValue,shocks,FORCS1(:,:,b),T,R,mv, mu);
end

mFORCS1 = mean(FORCS1,3);

tt = (1-ci)/2;
t1 = round(B*tt);
t2 = round(B*(1-tt));

forecasts.controled_variables = cV;
forecasts.instruments = cS;

for i = 1:EndoSize
    eval(['forecasts.cond.mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS1(i,:)'';']);
    tmp = sort(squeeze(FORCS1(i,:,:))');
    eval(['forecasts.cond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end

clear FORCS1;

FORCS2 = zeros(NumberOfStates,H+1,B);
for b=1:B
    FORCS2(:,1,b) = InitState;
end

randn('state',0);

for b=1:B
    shocks = sQ*randn(ExoSize,H);
    shocks(jdx,:) = zeros(length(jdx),H);
    FORCS2(:,:,b) = mcforecast3(0,H,mcValue,shocks,FORCS2(:,:,b),T,R,mv, mu);
end

mFORCS2 = mean(FORCS2,3);


for i = 1:EndoSize
    eval(['forecasts.uncond.mean.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ' = mFORCS2(i,:)'';']);
    tmp = sort(squeeze(FORCS2(i,:,:))');
    eval(['forecasts.uncond.ci.' deblank(M_.endo_names(oo_.dr.order_var(i),:)) ...
          ' = [tmp(t1,:)'' ,tmp(t2,:)'' ]'';']);
end

save('conditional_forecasts.mat','forecasts');