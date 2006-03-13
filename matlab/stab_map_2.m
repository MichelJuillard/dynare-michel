%function stab_map_2(x,alpha2,istab,fnam)
function stab_map_2(x,alpha2,fnam, ishock)
% function stab_map_2(x,alpha2,fnam, ishock)
%
% Copyright (C) 2005 Marco Ratto
% THIS PROGRAM WAS WRITTEN FOR MATLAB BY
% Marco Ratto,
% Unit of Econometrics and Statistics AF
% (http://www.jrc.cec.eu.int/uasa/),
% IPSC, Joint Research Centre
% The European Commission,
% TP 361, 21020 ISPRA(VA), ITALY
% marco.ratto@jrc.it 
%

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

if nargin<4,
    ishock=0;
end
ys_ = oo_.dr.ys;
dr_ = oo_.dr;
fname_ = M_.fname;
nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

c0=corrcoef(x);
c00=tril(c0,-1);
fig_nam_=[fname_,'_',fnam,'_corr_'];

ifig=0;
j2=0;
if ishock==0
    npar=estim_params_.np;
else
    npar=estim_params_.np+nshock;
end
for j=1:npar,
    i2=find(abs(c00(:,j))>alpha2);
    if length(i2)>0,
        for jx=1:length(i2),
            j2=j2+1;
            if mod(j2,12)==1,
                ifig=ifig+1;
                figure('name',['Correlations in the ',fnam,' sample ', num2str(ifig)]),
            end
            subplot(3,4,j2-(ifig-1)*12)
            %             bar(c0(i2,j)), 
            %             set(gca,'xticklabel',bayestopt_.name(i2)), 
            %             set(gca,'xtick',[1:length(i2)])
            %plot(stock_par(ixx(nfilt+1:end,i),j),stock_par(ixx(nfilt+1:end,i),i2(jx)),'.k')
            %hold on, 
            plot(x(:,j),x(:,i2(jx)),'.')
            %             xlabel(deblank(estim_params_.param_names(j,:)),'interpreter','none'), 
            %             ylabel(deblank(estim_params_.param_names(i2(jx),:)),'interpreter','none'), 
            if ishock,
                xlabel(bayestopt_.name{j},'interpreter','none'), 
                ylabel(bayestopt_.name{i2(jx)},'interpreter','none'), 
            else
                xlabel(bayestopt_.name{j+nshock},'interpreter','none'), 
                ylabel(bayestopt_.name{i2(jx)+nshock},'interpreter','none'), 
            end
            title(['cc = ',num2str(c0(i2(jx),j))])
            if (mod(j2,12)==0) & j2>0,
                saveas(gcf,[fig_nam_,int2str(ifig)])
            end
        end
    end
    if (j==(npar)) & j2>0,
        saveas(gcf,[fig_nam_,int2str(ifig)])
    end
    
end
if ifig==0,
    disp(['No correlation term >', num2str(alpha2),' found for ',fnam])
end
%close all
