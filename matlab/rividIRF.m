function rivid_irf(m, nroot, var_list_, iload)
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

%load caz
%M_.Sigma_e=caz;
if nargin<4, iload=0; end,
%%
MaxNumberOfPlotPerFigure = 9;% The square root must be an integer!
nn = sqrt(MaxNumberOfPlotPerFigure);
DirectoryName = CheckPath('Output');
MhDirectoryName = CheckPath('GSA');
MAX_nirfs = ceil(options_.MaxNumberOfBytes/(options_.irf*length(oo_.steady_state)*M_.exo_nbr)/8)+50;
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

load([ MhDirectoryName '/'  M_.fname '_prior'],'lpmat','istable')
x=lpmat(istable,:);
clear lpmat istable
NumberOfDraws=size(x,1);
B=NumberOfDraws; options_.B = B;

if ~iload
  try delete([MhDirectoryName '\' M_.fname '_IRFs*']);
  catch disp('No _IRFs files to be deleted!')
  end
  irun = 0;
  irun2 = 0;
  NumberOfIRFfiles = 1;
  ifil2 = 1;
  h = waitbar(0,'GSA (prior) IRFs...');
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
  u=[zeros(10,1); 1; zeros(options_.irf-1,1)];
  i0=[];
  izp=[];

  for b=1:B,
    zp=0;
    irun = irun+1;
    irun2 = irun2+1;
    deep = x(b,:);
    stock_param(irun2,:) = deep;
    set_parameters(deep);
    %set_coefficients(deep);  % only runs coefficients fixing the st. error of shocks
    dr = resol(oo_.steady_state,0);
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord) = M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    cs = transpose(chol(SS));
    for i = 1:M_.exo_nbr
      if SS(i,i) > 1e-13
        y=irf(dr,SS(M_.exo_names_orig_ord,i), options_.irf, options_.drop,options_.replic,options_.order);
        if options_.relative_irf
          y = 100*y/cs(i,i);
        end
        for j = 1:M_.endo_nbr,
          jj=strmatch(deblank(M_.endo_names(j,:)),var_list_,'exact');
          if jj,
            yy=[zeros(10,1); y(j,:)'];
            tf=squeeze(m(i,jj,:))';
            [th,stats]=riv([yy u],tf,10);
            if stats(3)<1-1.e-6,
              disp('');
              %[th,stats]=rivid([yy u],[[tf(1:2)-1 tf(3:4)]; [tf(1:2)+1 tf(3:4)]],1,10);
              [a,bb,c,P,d]=getpar(th);
              tf(1)=length(a)-1;
              tf(3)=length(find(bb==0));
              tf(2)=length(find(bb));
              tf(4)=0;
            end
            sta(i,jj,b,:)=stats;
            [a,bb,c,P,d]=getpar(th);
            am{i,jj}(b,1:tf(1))=a(2:end);
            bm{i,jj}(b,1:tf(2))=bb(tf(3)+1:end);
            ar=roots(a);
            for ir=1:tf(1),
              if find(abs(roots(bb)-ar(ir))<1.e-5)
                zp=1;
                [th,stats]=riv([yy u],[tf(1)-1,tf(2)-1,tf(3:end)],10);
                sta(i,jj,b,:)=stats;
                [a,bb,c,P,d]=getpar(th);
                am{i,jj}(b,:)=am{i,jj}(b,:)*0;
                bm{i,jj}(b,:)=bm{i,jj}(b,:)*0;
                am{i,jj}(b,1:tf(1)-1)=a(2:end);
                bm{i,jj}(b,1:tf(2)-1)=bb(tf(3)+1:end);
              end
            end
          end
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
    if irun2 == MAX_nruns | b == B
      if b == B
        stock_param = stock_param(1:irun2,:);
      end
      stock = stock_param;
      save([MhDirectoryName '/' M_.fname '_param_irf' int2str(ifil2)],'stock');
      ifil2 = ifil2 + 1;
      irun2 = 0;
    end
    if min(min(sta(:,:,b,3)))<(1-1.e-6), % | zp,
      i0=[i0 b];
      if zp,
        izp=[izp b];
      end
    end,

    waitbar(b/B,h);
  end
  NumberOfIRFfiles = NumberOfIRFfiles-1;
  ifil2 = ifil2-1;
  close(h);
  ReshapeMatFiles('irf','gsa');
  i2=find(~ismember([1:B],i0));
  save([MhDirectoryName '/' M_.fname '_rivid'],'x','am','bm','i0','i2','izp','sta');
else
  load([MhDirectoryName '/' M_.fname '_rivid'],'x','am','bm','i0','i2','izp','sta');
end

a0=[];
for i=1:M_.exo_nbr,
  for j=1:size(var_list_,1)
    a0=[a0 am{i,j}(i2,:)];
  end
end


r0={};
r00=[];
a00=[];
%nroot=3;
izp=[];
vzp={};
h = waitbar(0,'Analysing roots of denominators...');

for b=1:length(i2),
  roo=[];
  zp=0;
  wzp=[];
  for i=1:M_.exo_nbr,
    for j=1:size(var_list_,1),
      tf=squeeze(m(i,j,:))';
      bb=bm{i,j}(i2(b),:);
      r0{i,j}(b,:)=sort(roots([1 am{i,j}(i2(b),:)]));
      %       bb=bm{i,j}(b,:);
      %       r0{i,j}(b,:)=sort(roots([1 am{i,j}(b,:)]));
      for ir=1:tf(1),
        if r0{i,j}(b,ir)~=0 & find(abs(roots(bb)-r0{i,j}(b,ir))<1.e-5)
          zp=1;
          wzp=[wzp; [i j]];
        end
      end

      roo=[roo r0{i,j}(b,:)];
      %       if j>1 | i>1,
      %         f0=[];
      %       for l=1:m(i,j,1)
      %         f0(l)=find(abs(((r0{i,j}(b,:))-(r0{i,1}(b,l))))<=1.e-6);
      %       end
      %       r0{i,j}(b,:)=r0{i,j}(b,f0);
      %       end
    end
  end
  roo=sort(roo);
  crit=1;
  rxx=[];
  while length(rxx)~=nroot & zp==0 & crit>1.e-10,
    crit=crit*0.5;
    ix=find(abs(diff(sign(real(roo)).*abs(roo)))>crit);
    ix=[1 ix+1];
    rx=roo(ix);
    rxx=[];
    for jx=1:length(rx),
      rxx=[rxx rx(jx)];
      if find(imag(rx(jx))),
        rxx=[rxx conj(rx(jx))];
      end
    end
  end
  if zp==0 & crit>1.e-10,
    r00=[r00; rxx];
    a00=[a00; poly(rxx)];
  else
    izp=[izp i2(b)];
    %     izp=[izp b];
    vzp=[vzp; {wzp}];
    r00=[r00; zeros(1,nroot)];
    a00=[a00; poly(zeros(1,nroot))];
  end
  waitbar(b/length(i2),h)
  %   waitbar(b/1669,h)
end
close(h)

i1=find(~ismember([i2],izp));
save([MhDirectoryName '/' M_.fname '_rivid'],'r00','a00','-append');


if options_.opt_gsa.morris,
  load([MhDirectoryName '/' M_.fname '_prior'],'lpmat');
  lmat=[];
  inx=[];
  for j=1:options_.opt_gsa.morris_ntra,
    z(j)=length(find(ismember([1+(estim_params_.np+1)*(j-1):(estim_params_.np+1)*j],izp)));
    if z(j)==0,
      inx=[inx [1+12*(j-1):12*j]];
      lmat=[lmat; lpmat([1+12*(j-1):12*j],:)];
    end
  end,
  nliv = options_.opt_gsa.morris_nliv;
  j0=0;
  for j=1:nroot,
    if (max(a00(inx,j+1))-min(a00(inx,j+1)))>1.e-10,
      j0=j0+1;
      [SAmeas, SAMorris(:,:,j0)] = Morris_Measure_Groups(estim_params_.np, lmat, a00(inx,j+1),nliv);
    end
  end
  j1=j0;
  for j=1:M_.exo_nbr,
    SAMo=[];
    j00=0;
    yb=[];
    y=am{j,1};
    for l=1:size(y,2),
      y0=y(:,l);
      if max(y0(inx))-min(y0(inx))>1.e-10,
        j00=j00+1;
        yb=[yb y0(inx)];
        [SAmeas, SAMo(:,:,j00)] = Morris_Measure_Groups(estim_params_.np, lmat, y0(inx),nliv);
      end
    end
    for i=1:size(var_list_,1),
      y=bm{j,i};
      for l=1:size(y,2),
        y0=y(:,l);
        if max(y0(inx))-min(y0(inx))>1.e-10,
          j00=j00+1;
          yb=[yb y0(inx)];
          [SAmeas, SAMo(:,:,j00)] = Morris_Measure_Groups(estim_params_.np, lmat, y0(inx),nliv);
          j0=j0+1;
          [SAmeas, SAMorris(:,:,j0)] = Morris_Measure_Groups(estim_params_.np, lmat, y0(inx),nliv);
        end
      end
    end
    SAM = squeeze(SAMo(:,1,:));
    SAno=[];
    for jj=1:j00
      SAno(:,jj)=SAM(:,jj)./max(SAM(:,jj));
    end
    ybb(j)={yb};
    figure, bar(SAno), title(M_.exo_names(j,:))
  end

  SAM = squeeze(SAMorris(:,1,:));
  for j=1:j0,
    SAnorm(:,j)=SAM(:,j)./max(SAM(:,j));
  end

  figure, bar(SAnorm), title('All')
else
  gsa_flag=0; %0; %-2 for GSA estimation
  x=x(i1,:);
  ya=a00(i1,2:end);
  j0=size(ya,2);
  for j=1:M_.exo_nbr,
    j00=0;
    yb=[];
    yr=[];
    y=am{j,1}(i1,:);
    for l=1:size(y,2),
      y0=y(:,l);
      if max(y0)-min(y0)>1.e-10,
        j00=j00+1;
        yb=[yb y0];
      end
    end
    for i=1:size(var_list_,1),
      y=bm{j,i}(i1,:);
      for l=1:size(y,2),
        y0=y(:,l);
        if max(y0)-min(y0)>1.e-10,
          j00=j00+1;
          yb=[yb y0];
          j0=j0+1;
          ya=[ya y0];
        end
      end
    end
    ybb(j)={yb};
    for jj=1:j00,
      [dum, is]=sort(yb(:,jj));
      yr(is,jj)=[1:length(i1)]'./length(i1);
    end
    ybr(j)={yr};
  end
  yr=[];
  for j=1:j0,
    [dum, is]=sort(ya(:,j));
    yr(is,j)=[1:length(i1)]'./length(i1);
  end
  RividDir = CheckPath('GSA/rivid');
  for j=1:j0,
  gsa_(j) = gsa_sdp_fn(ya(:,j), x, gsa_flag, [RividDir,'/map',int2str(j)], M_.param_names);
  end
%   for j=1:j0,
%   gsa_y(j) = gsa_sdp_fn(ya(:,j), ya(:,[1:j-1,j+1:j0]), 0, [RividDir,'/map_y',int2str(j)], M_.param_names);
%   end
  for j=1:j0, S(j,:)=gsa_(j).multivariate.si; end,
  [v,d]=eig(corrcoef(ya));
  ee=diag(d);
  %id=find((er./length(er))>0.01);
  id=find(cumsum(ee(end:-1:1))./j0<0.99);
  id=length(id)+1;
  jpc=0;
  for j=j0:-1:(j0-id(1)+1),
    jpc=jpc+1;
    gsa_pca(jpc) = gsa_sdp_fn(ya*v(:,j), x, gsa_flag, [RividDir,'/map_pca',int2str(jpc)], M_.param_names);
  end
  for j=1:jpc, S_pca(j,:)=gsa_pca(j).multivariate.si; end,
  for j=1:j0,
  gsar_(j) = gsa_sdp_fn(yr(:,j), x, gsa_flag, [RividDir,'/map_rank',int2str(j)], M_.param_names);
  end
  for j=1:j0, Sr(j,:)=gsar_(j).multivariate.si; end,
  [v,d]=eig(corrcoef(yr));
  er=diag(d);
  %id=find((er./length(er))>0.01);
  id=find(cumsum(er(end:-1:1))./j0<0.99);
  id=length(id)+1;
  jpc=0;
  for j=j0:-1:(j0-id(1)+1),
    jpc=jpc+1;
    gsar_pca(jpc) = gsa_sdp_fn(yr*v(:,j), x, gsa_flag, [RividDir,'/map_rank_pca',int2str(jpc)], M_.param_names);
  end
  for j=1:jpc, Sr_pca(j,:)=gsar_pca(j).multivariate.si; end,
  fi=zeros(id,estim_params_.np);
  for j=1:id,
    yup=find(Sr_pca(j,:)>0.05);
    if ~isempty(yup),
      fi(j,yup)=ones(1,length(yup));
    end
  end
%   for j=1:j0-id,
%     gsar_pca_0(j) = gsa_sdp_fn(yr*v(:,j), x, gsa_flag, [RividDir,'/rivid/map_rank_pca_0',int2str(j)], M_.param_names);
%   end
%   for j=1:j0-id, Sr_pca_0(j,:)=gsar_pca_0(j).multivariate.si; end,
%   yr0=zeros(length(i1),1);
%   for j=1:10, %j0-id,
%     yr0 = yr0+yr*v(:,j);
%   end
%   gsar_pca_00 = gsa_sdp_fn(yr0, x, gsa_flag, [RividDir,'/map_rank_pca_00'], M_.param_names);
%   Sr_pca_00=gsar_pca_00.multivariate.si;
%   y0 = getIRFRuns(1,1,'gsa');
%   prc=[2.5 16 84 97.5];
%   [paths_avg,res,bounds,res1]=probam(y0',prc,2);
%   figure,
%   plot([1:options_.irf],paths_avg,'k-',[1:options_.irf],bounds(1,:),'k:', [1:options_.irf], bounds(2,:),'k:');
  save([MhDirectoryName '/' M_.fname '_rivid'],'v','d','-append');
jpc=0;
for j=j0:-1:(j0-id(1)+1),
  jpc=jpc+1;
  y0=yr*v(:,j);
  rbf_1(jpc)=rbf(x, y0, 1);
  disp(['Output ',num2str(jpc),' out of ',num2str(id(1)),' done!'])

end
  jpc=0;
  
  for j=j0:-1:(j0-id(1)+1),
    jpc=jpc+1;
    y0=yr*v(:,j);
    ivec=find(gsar_pca(jpc).multivariate.si>0.01);
    y0=y0-sum(gsar_pca(jpc).multivariate.f(:,ivec),2);
    rbf_2(jpc)=rbf(x(:,ivec), y0, 2);
    disp(['Output ',num2str(jpc),' out of ',num2str(id(1)),' done!'])

  end
  save([MhDirectoryName '/' M_.fname '_rivid'],'rbf_1','rbf_2','-append');

end
