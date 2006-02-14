function stab_map_2(x,alpha2,istab)
% function stab_map_2(x,alpha2,istab)
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

global bayestopt_ estim_params_ dr_ options_ ys_ fname_

c0=corrcoef(x);
c00=tril(c0,-1);
if istab,
    fig_nam_=[fname_,'_stab_corr_'];
else
    fig_nam_=[fname_,'_unstab_corr_'];
end

ifig=0;
j2=0;
npar=estim_params_.np;
for j=1:npar,
    i2=find(abs(c00(:,j))>alpha2);
    if length(i2)>0,
        for jx=1:length(i2),
            j2=j2+1;
            if mod(j2,12)==1,
                ifig=ifig+1;
                if istab
                    figure('name',['Correlations in the stable sample ', num2str(ifig)]),
                else
                    figure('name',['Correlations in the unstable sample ', num2str(ifig)]),
                end
            end
            subplot(3,4,j2-(ifig-1)*12)
            %             bar(c0(i2,j)), 
            %             set(gca,'xticklabel',bayestopt_.name(i2)), 
            %             set(gca,'xtick',[1:length(i2)])
            %plot(stock_par(ixx(nfilt+1:end,i),j),stock_par(ixx(nfilt+1:end,i),i2(jx)),'.k')
            %hold on, 
            plot(x(:,j),x(:,i2(jx)),'.')
            xlabel(deblank(estim_params_.param_names(j,:)),'interpreter','none'), 
            ylabel(deblank(estim_params_.param_names(i2(jx),:)),'interpreter','none'), 
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
%close all
