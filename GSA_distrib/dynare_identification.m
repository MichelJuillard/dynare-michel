function [pdraws, idemodel, idemoments] = dynare_identification()

% main 

global M_ options_ oo_ bayestopt_ estim_params_


options_ = set_default_option(options_,'datafile',[]);
options_.mode_compute = 0;
[data,rawdata]=dynare_estimation_init([],1);
% computes a first linear solution to set up various variables
dynare_resolve;

options_.prior_mc=2000;

SampleSize = options_.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

prior_draw(1,bayestopt_);
IdentifDirectoryName = CheckPath('identification');

indx = estim_params_.param_vals(:,1);
indexo=[];
if ~isempty(estim_params_.var_exo)
  indexo = estim_params_.var_exo(:,1);
end
useautocorr = 0;
nlags = 3;

iteration = 0;
loop_indx = 0;

h = waitbar(0,'Monte Carlo identification checks ...');

while iteration < SampleSize,
  loop_indx = loop_indx+1;
  params = prior_draw();
  set_all_parameters(params);

  [JJ, H] = getJJ(M_,oo_,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);  
  
  if ~isempty(JJ),
    iteration = iteration + 1;
    pdraws(iteration,:) = params';
    [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), ...
      idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), ...
      idemodel.cond(iteration), idemoments.cond(iteration), ...
      idemodel.ee(:,iteration), idemoments.ee(:,iteration), ...
      idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
      idemodel.indno{iteration}, idemoments.indno{iteration}] = ...
      identification_checks(H,JJ, bayestopt_);   
    
      waitbar(iteration/SampleSize,h)
  end
end

close(h)

npar = size(pdraws,2);
jok = 0;
jokP = 0;
jokJ = 0;
jokPJ = 0;
for j=1:npar,
  if any(idemodel.ind(j,:)==0),
    pno = 100*length(find(idemodel.ind(j,:)==0))/SampleSize;
    disp(['Parameter ',bayestopt_.name{j},' is not identified in the model for ',num2str(pno),'% of MC runs!' ])
    disp(' ')
  end
  if any(idemoments.ind(j,:)==0),
    pno = 100*length(find(idemoments.ind(j,:)==0))/SampleSize;
    disp(['Parameter ',bayestopt_.name{j},' is not identified by J moments for ',num2str(pno),'% of MC runs!' ])
    disp(' ')
  end
  if any(idemodel.ind(j,:)==1),
    iok = find(idemodel.ind(j,:)==1);
    jok = jok+1;
    kok(jok) = j;
    mmin(jok,1) = min(idemodel.Mco(j,iok));
    mmean(jok,1) = mean(idemodel.Mco(j,iok));
    mmax(jok,1) = max(idemodel.Mco(j,iok));
    [ipmax, jpmax] = find(abs(squeeze(idemodel.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
    if ~isempty(ipmax)
    jokP = jokP+1;
    kokP(jokP) = j;
    ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
    [N,X]=hist(ipmax,[1:npar]);
    jpM(jokP)={find(N)};
    NPM(jokP)={N(find(N))./SampleSize.*100};
    pmeanM(jokP)={mean(squeeze(idemodel.Pco(j,find(N),iok))')};
    pminM(jokP)={min(squeeze(idemodel.Pco(j,find(N),iok))')};
    pmaxM(jokP)={max(squeeze(idemodel.Pco(j,find(N),iok))')};
    end
  end
  if any(idemoments.ind(j,:)==1),
    iok = find(idemoments.ind(j,:)==1);
    jokJ = jokJ+1;
    kokJ(jokJ) = j;
    mminJ(jokJ,1) = min(idemoments.Mco(j,iok));
    mmeanJ(jokJ,1) = mean(idemoments.Mco(j,iok));
    mmaxJ(jokJ,1) = max(idemoments.Mco(j,iok));
    [ipmax, jpmax] = find(abs(squeeze(idemoments.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
    if ~isempty(ipmax)
    jokPJ = jokPJ+1;
    kokPJ(jokPJ) = j;
    ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
    [N,X]=hist(ipmax,[1:npar]);
    jpJ(jokPJ)={find(N)};
    NPJ(jokPJ)={N(find(N))./SampleSize.*100};
    pmeanJ(jokPJ)={mean(squeeze(idemoments.Pco(j,find(N),iok))')};
    pminJ(jokPJ)={min(squeeze(idemoments.Pco(j,find(N),iok))')};
    pmaxJ(jokPJ)={max(squeeze(idemoments.Pco(j,find(N),iok))')};
    end
  end
end

dyntable('Multi collinearity in the model:',strvcat('param','min','mean','max'), ...
  strvcat(bayestopt_.name(kok)),[mmin, mmean, mmax],10,10,6);

dyntable('Multi collinearity for moments in J:',strvcat('param','min','mean','max'), ...
  strvcat(bayestopt_.name(kokJ)),[mminJ, mmeanJ, mmaxJ],10,10,6);

for j=1:length(kokP),
dyntable([bayestopt_.name{kokP(j)},' pairwise correlations in the model'],strvcat(' ','min','mean','max'), ...
  strvcat(bayestopt_.name{jpM{j}}),[pminM{j}' pmeanM{j}' pmaxM{j}'],10,10,3);  
end

for j=1:length(kokPJ),
dyntable([bayestopt_.name{kokPJ(j)},' pairwise correlations in J moments'],strvcat(' ','min','mean','max'), ...
  strvcat(bayestopt_.name{jpJ{j}}),[pminJ{j}' pmeanJ{j}' pmaxJ{j}'],10,10,3);  
end

save([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments')



