function redform_map(anamendo, anamlagendo, anamexo, iload, pprior, ilog, threshold, ksstat, alpha2, dirname)
%function redform_map(namendo, namlagendo, namexo, iload, pprior, ilog, threshold, ksstat, alpha2, dirname)
% copyright Marco Ratto 2006
global M_ oo_ estim_params_ options_

pnames = M_.param_names(estim_params_.param_vals(:,1),:);
if nargin<4 | isempty(iload),
  iload=0;
end
if nargin<5 | isempty(pprior),
  pprior=1;
end
if nargin<6 | isempty(ilog),
  ilog=0;
end
if nargin<7,
  threshold=[];
end
if nargin<8,
  ksstat=0.1;
end
if nargin<9,
  alpha2=0.2;
end
if nargin<10,
  dirname='';
end

if pprior
  load([dirname,'\',M_.fname,'_prior']);
  adir=[dirname '\redform_stab'];
else
  load([dirname,'\',M_.fname,'_mc']);
  adir=[dirname '\redform_mc'];
end
if isempty(dir(adir))
  mkdir(adir)
end
adir0=pwd;
%cd(adir)

nspred=size(T,2)-M_.exo_nbr;
x0=lpmat(istable,:);
[kn, np]=size(x0);
nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));
for j=1:size(anamendo,1)
  namendo=deblank(anamendo(j,:));
  iendo=strmatch(namendo,M_.endo_names(oo_.dr.order_var,:),'exact');
  for jx=1:size(anamexo,1)
    namexo=deblank(anamexo(jx,:));
    iexo=strmatch(namexo,M_.exo_names,'exact');
    
    if ~isempty(iexo),
      %y0=squeeze(T(iendo,iexo+nspred,istable));
      y0=squeeze(T(iendo,iexo+nspred,:));
      clear lpmat lpmat0 egg iunstable yys
      xdir0 = [adir,'\',namendo,'_vs_', namexo];
      if isempty(threshold)
        redform_private(x0, y0, iload, pnames, namendo, namexo, xdir0)
      else
        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
        xdir = [xdir0,'_cut'];
        redform_private(x0(iy,:), y0(iy), iload, pnames, namendo, namexo, xdir)
        delete([xdir, '\*cut*.*'])
        [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
        indsmirnov = find(dproba>ksstat);
        stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
        stab_map_2(x0(iy,:),alpha2,'cut',xdir)
        stab_map_2(x0(iyc,:),alpha2,'trim',xdir)
      end
      if ilog==-1,
        yy=log(-y0);
        xdir=[xdir0,'_minuslog'];
        redform_private(x0, yy, iload, pnames, namendo, namexo, xdir)
      elseif ilog==1,
        yy=log(y0);
        xdir=[xdir0,'_log'];
        redform_private(x0, yy, iload, pnames, namendo, namexo, xdir)
      end
    end
  end
  for je=1:size(anamlagendo,1)
    namlagendo=deblank(anamlagendo(je,:));
    ilagendo=strmatch(namlagendo,M_.endo_names(oo_.dr.order_var(oo_.dr.nstatic+1:oo_.dr.nstatic+nsok),:),'exact');
    
    if ~isempty(ilagendo),
      %y0=squeeze(T(iendo,ilagendo,istable));
      y0=squeeze(T(iendo,ilagendo,:));
      xdir0 = [adir,'\',namendo,'_vs_', namlagendo];
      if isempty(threshold)
        redform_private(x0, y0, iload, pnames, namendo, namlagendo, xdir0)
      else
        iy=find( (y0>threshold(1)) & (y0<threshold(2)));
        iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
        xdir = [xdir0,'_cut'];
        redform_private(x0(iy,:), y0(iy), iload, pnames, namendo, namlagendo, xdir)
        delete([xdir, '\*cut*.*'])
        [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
        indsmirnov = find(dproba>ksstat);
        stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
        stab_map_2(x0(iy,:),alpha2,'cut',xdir)
        stab_map_2(x0(iyc,:),alpha2,'trim',xdir)
      end
      if ilog==-1,
        yy=log(-y0);
        xdir=[xdir0,'_minuslog'];
        redform_private(x0, yy, iload, pnames, namendo, namlagendo, xdir)
      elseif ilog==1,
        yy=log(y0);
        xdir=[xdir0,'_log'];
        redform_private(x0, yy, iload, pnames, namendo, namlagendo, xdir)
      end
    end
  end
end


function redform_private(x0, y0, iload, pnames, namy, namx, xdir)
figure, hist(y0,30), title([namy,' vs. ', namx])
if isempty(dir(xdir))
  mkdir(xdir)
end
saveas(gcf,[xdir,'\', namy,'_vs_', namx])
eval(['print -depsc2 ' xdir,'\', namy,'_vs_', namx]);
eval(['print -dpdf ' xdir,'\', namy,'_vs_', namx]);
fname=[xdir,'\map'];
if iload==0,
  gsa_ = gsa_sdp_fn(y0, x0, -2, fname, pnames, 1);
else
  gsa_ = gsa_sdp_fn(y0, x0, 0, fname, pnames, 1);
end
