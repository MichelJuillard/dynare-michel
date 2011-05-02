function plot_identification(params,idemoments,idehess,idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName, save_figure)
% function plot_identification(params,idemoments,idehess,idemodel, idelre, advanced, tittxt, name, IdentifDirectoryName, save_figure)

% Copyright (C) 2008-2011 Dynare Team
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

global M_ options_

if nargin<10 || isempty(save_figure),
    save_figure=0;
end

[SampleSize, nparam]=size(params);
siJnorm = idemoments.siJnorm;
siHnorm = idemodel.siHnorm;
siLREnorm = idelre.siLREnorm;

% if prior_exist,
%     tittxt = 'Prior mean - ';
% else
%     tittxt = '';
% end

if SampleSize == 1,
    siJ = idemoments.siJ;
    normJ = max(abs(siJ)')';
    figure('Name',[tittxt, 'Identification using info from observables']),
    subplot(211)
    mmm = (idehess.ide_strength_J);
    [ss, is] = sort(mmm);
    bar(log([idehess.ide_strength_J(:,is)' idehess.ide_strength_J_prior(:,is)']))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    legend('relative to param value','relative to prior std','Location','Best')
    if  idehess.flag_score,
        title('Identification strength in the asymptotic Information matrix (log-scale)')
    else
        title('Identification strength in the moments (log-scale)')
    end
    
    subplot(212)
    
    mmm = (siJnorm)'./max(siJnorm);
    if advanced,
        mmm1 = (siHnorm)'./max(siHnorm);
        mmm=[mmm mmm1];
        mmm1 = (siLREnorm)'./max(siLREnorm);
        offset=length(siHnorm)-length(siLREnorm);
        mmm1 = [NaN(offset,1); mmm1];
        mmm=[mmm mmm1];
    end        
        
    bar(mmm(is,:))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    if advanced,
        legend('Moments','Model','LRE model','Location','Best')
    end
    title('Sensitivity bars')
    
    if advanced
        % identificaton patterns
        for  j=1:size(idemoments.cosnJ,2),
            pax=NaN(nparam,nparam);
            fprintf('\n\n')
            disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
            fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
            for i=1:nparam,
                namx='';
                for in=1:j,
                    dumpindx = idemoments.pars{i,j}(in);
                    if isnan(dumpindx),
                        namx=[namx ' ' sprintf('%-15s','--')];
                    else
                        namx=[namx ' ' sprintf('%-15s',name{dumpindx})];
                        pax(i,dumpindx)=idemoments.cosnJ(i,j);
                    end
                end
                fprintf('%-15s [%s] %10.3f\n',name{i},namx,idemoments.cosnJ(i,j))
            end
            figure('name',[tittxt,'Collinearity patterns with ', int2str(j) ,' parameter(s)']),
            imagesc(pax,[0 1]);
            set(gca,'xticklabel','')
            set(gca,'yticklabel','')
            for ip=1:nparam,
                text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
                text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
            end
            colorbar;
            ax=colormap;
            ax(1,:)=[0.9 0.9 0.9];
            colormap(ax);
            if nparam>10,
                set(gca,'xtick',(5:5:nparam))
                set(gca,'ytick',(5:5:nparam))
            end
            set(gca,'xgrid','on')
            set(gca,'ygrid','on')
            if save_figure
                saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_collinearity_', int2str(j)])
                eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_collinearity_', int2str(j)]);
                eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_collinearity_', int2str(j)]);
                if options_.nograph, close(gcf); end
            end
        end
        disp('')
        if idehess.flag_score,
            [U,S,V]=svd(idehess.AHess,0);
            if nparam<5,
                f1 = figure('name',[tittxt,'Identification patterns (Information matrix)']);
            else
                f1 = figure('name',[tittxt,'Identification patterns (Information matrix): SMALLEST SV']);
                f2 = figure('name',[tittxt,'Identification patterns (Information matrix): HIGHEST SV']);
            end
        else
            [U,S,V]=svd(siJ./normJ(:,ones(nparam,1)),0);
            if nparam<5,
                f1 = figure('name',[tittxt,'Identification patterns (moments)']);
            else
                f1 = figure('name',[tittxt,'Identification patterns (moments): SMALLEST SV']);
                f2 = figure('name',[tittxt,'Identification patterns (moments): HIGHEST SV']);
            end
        end
        for j=1:min(nparam,8),
            if j<5,
                figure(f1),
                jj=j;
            else
                figure(f2),
                jj=j-4;
            end
            subplot(4,1,jj),
            if j<5
                bar(abs(V(:,end-j+1))),
                Stit = S(end-j+1,end-j+1);
            else
                bar(abs(V(:,jj))),
                Stit = S(jj,jj);
            end
            set(gca,'xticklabel','')
            if j==4 || j==nparam || j==8,
                for ip=1:nparam,
                    text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
                end
            end
            title(['Singular value ',num2str(Stit)])
        end
        if save_figure,
            figure(f1);
            saveas(f1,[IdentifDirectoryName,'/',M_.fname,'_ident_pattern_1'])
            eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_pattern_1']);
            eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_pattern_1']);
            if nparam>4,
                figure(f2),
                saveas(f2,[IdentifDirectoryName,'/',M_.fname,'_ident_pattern_2'])
                eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_pattern_2']);
                eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_pattern_2']);
            end
        end
    end
    
else
    figure('Name',['MC sensitivities']),
    mmm = (idehess.ide_strength_J);
    [ss, is] = sort(mmm);
    mmm = mean(siJnorm)';
    mmm = mmm./max(mmm);
    if advanced,
        mmm1 = mean(siHnorm)';
        mmm=[mmm mmm1./max(mmm1)];
        mmm1 = mean(siLREnorm)';
        offset=size(siHnorm,2)-size(siLREnorm,2);
        mmm1 = [NaN(offset,1); mmm1./max(mmm1)];
        mmm=[mmm mmm1];
    end        
        
    bar(mmm(is,:))
    set(gca,'xlim',[0 nparam+1])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    if advanced,
        legend('Moments','Model','LRE model','Location','Best')
    end
    title('MC mean of sensitivity measures')
    if advanced,
        options_.nograph=1;
        figure('Name','MC Condition Number'),
        subplot(221)
        hist(log10(idemodel.cond))
        title('log10 of Condition number in the model')
        subplot(222)
        hist(log10(idemoments.cond))
        title('log10 of Condition number in the moments')
        subplot(223)
        hist(log10(idelre.cond))
        title('log10 of Condition number in the LRE model')
        saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_COND'])
        eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_COND']);
        eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_COND']);
        if options_.nograph, close(gcf); end
        ncut=floor(SampleSize/10*9);
        [~,is]=sort(idelre.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberLRE', 1, [], IdentifDirectoryName, 0.1);
        [~,is]=sort(idemodel.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberModel', 1, [], IdentifDirectoryName, 0.1);
        [~,is]=sort(idemoments.cond);
        [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), 'MC_HighestCondNumberMoments', 1, [], IdentifDirectoryName, 0.1);
%         [proba, dproba] = stab_map_1(idemoments.Mco', is(1:ncut), is(ncut+1:end), 'HighestCondNumberMoments_vs_Mco', 1, [], IdentifDirectoryName);
        for j=1:nparam,
%             ibeh=find(idemoments.Mco(j,:)<0.9);
%             inonbeh=find(idemoments.Mco(j,:)>=0.9);
%             if ~isempty(ibeh) && ~isempty(inonbeh)
%                 [proba, dproba] = stab_map_1(params, ibeh, inonbeh, ['HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName);
%             end
            [~,is]=sort(idemoments.Mco(:,j));
            [proba, dproba] = stab_map_1(params, is(1:ncut), is(ncut+1:end), ['MC_HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName, 0.15);
        end
    end
end

% disp_identification(params, idemodel, idemoments, name)