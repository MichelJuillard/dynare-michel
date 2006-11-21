function plot_IRF(xname, vargin, s1, texname)
global M_ oo_ options_

fname_ =M_.fname;
DirectoryName = CheckPath('Output');

if nargin<4,
  for j=1:M_.endo_nbr,
    texname{j}=deblank(M_.endo_names_tex(j,:));
  end
end
  

iendo=[];
for i=1:length(vargin)
  iendo=[iendo, strmatch(vargin{i},M_.endo_names,'exact')];
end
iexo=[];
for i=1:length(xname)
  iexo=[iexo, strmatch(xname{i},M_.exo_names,'exact')];
end

for i=1:length(iexo),
  nfig=0;
  nplo=0;
  for j=1:length(vargin),
    
    name = [vargin{j} '_' deblank(M_.exo_names(iexo(i),:))];
    eval(['MeanIRF=s1.' name,';']);

    if max(abs(MeanIRF)) > 1e-6 ,
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))])
        nfig=nfig+1;
      end
      subplot(3,3,nplo)
      plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
      hold on,
      plot(1:options_.irf,MeanIRF(:),'-k','linewidth',1)
      xlim([1 options_.irf]);
      hold off
      if options_.TeX,
        if nargin<4,
          title(texname{iendo(j)},'interpreter','tex')
        else
          title(texname{j},'interpreter','tex')
        end          
      else
        title(vargin{j},'interpreter','none')
      end
    end
    if (mod(nplo,9)==0 | j==length(vargin)) & nplo,
      saveas(gcf,[DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName,'\',fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
      close(gcf)
      nplo=0;
    end
  end
end
