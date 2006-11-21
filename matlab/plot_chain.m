function plot_chain
global M_ options_ estim_params_ bayestopt_

npar = estim_params_.np+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nvx;

nblck = options_.mh_nblck; 
DirectoryName = CheckPath('metropolis');
OutDirectoryName = CheckPath('output');
load([DirectoryName '/'  M_.fname '_mh_history'])
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));

ifig=0;
for j=1:npar,
  if mod(j,9)==1,
    ifig=ifig+1;
    iplot=0;
    hh=figure('Name',['Markov Chain plot ',int2str(ifig)]);
  end
  iplot=iplot+1;
  Draws = GetAllPosteriorDraws(j,1,1,TotalNumberOfMhFiles,TotalNumberOfMhDraws);
  Draws=reshape(Draws,TotalNumberOfMhDraws,nblck);
  Draws = stand_(Draws);
  subplot(3,3,iplot)
  plot(Draws)
  hold on, plot(cumsum(Draws)./[1:length(Draws)]','r')
  title(bayestopt_.name{j},'interpreter','none')
  if mod(j,9)==0 | j==npar,
    saveas(hh,[OutDirectoryName '\' M_.fname '_chain_' int2str(ifig)])
    eval(['print -depsc2 ' OutDirectoryName '/'  M_.fname '_chain_' int2str(ifig)]);
    eval(['print -dpdf ' OutDirectoryName '/' M_.fname  '_chain_' int2str(ifig)]);
    close(hh)
  end
  
end


