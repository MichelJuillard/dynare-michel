function McMCDiagnostics(Origin,alpha,NumberOfParameters,NumberOfSimulations,LastFileNumber,NumberOfBlocks,...
			 NumberOfSimulationsPerFile,CumulatedNumberOfSimulationsPerFile)
  
global M_ options_

StepSize = ceil((NumberOfSimulations-Origin)/100); 
time = 1:NumberOfSimulations;
xx = Origin:StepSize:NumberOfSimulations;
NumberOfLines = length(xx);
tmp = zeros(NumberOfSimulations*NumberOfBlocks,3);
UDIAG = zeros(NumberOfLines,6,NumberOfParameters);
if NumberOfSimulations < Origin
  error('MH: The number of simulations is to small to compute the MCMC convergence diagnostics.')
end

disp('MH: Univariate convergence diagnostic, Brooks and Gelman (1998):')
for j=1:NumberOfParameters
  fprintf('    Parameter %d...  ',j);
  for b = 1:NumberOfBlocks
    startline = 0;
    for file = 0:LastFileNumber
      eval(['load ' M_.fname '_mh' int2str(file) '_blck' int2str(b)]);
      clear logpo2 post2;
      tmp((b-1)*NumberOfSimulations+startline+1:...
	  (b-1)*NumberOfSimulations+CumulatedNumberOfSimulationsPerFile(file+1),1) = x2(:,j);
      clear x2;
      startline = startline+NumberOfSimulationsPerFile(file+1);
    end	
  end
  tmp(:,2) = kron((1:NumberOfBlocks)',ones(NumberOfSimulations,1));
  tmp(:,3) = kron(ones(NumberOfBlocks,1),time'); 
  tmp = sortrows(tmp,1);
  ligne   = 0;
  for iter  = Origin:StepSize:NumberOfSimulations
    ligne = ligne+1;
    linea = ceil(0.5*iter);
    n     = iter-linea+1;
    cinf  = round(n*alpha/2);
    csup  = round(n*(1-alpha/2));
    CINF  = round(NumberOfBlocks*n*alpha/2);
    CSUP  = round(NumberOfBlocks*n*(1-alpha/2));
    temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
    UDIAG(ligne,1,j) = temp(CSUP,1)-temp(CINF,1);
    moyenne = mean(temp(:,1));%% Pooled mean.
    UDIAG(ligne,3,j) = sum((temp(:,1)-moyenne).^2)/(NumberOfBlocks*n-1);
    UDIAG(ligne,5,j) = sum(abs(temp(:,1)-moyenne).^3)/(NumberOfBlocks*n-1);
    for i=1:NumberOfBlocks
      pmet = temp(find(temp(:,2)==i));
      UDIAG(ligne,2,j) = UDIAG(ligne,2,j) + pmet(csup,1)-pmet(cinf,1);
      moyenne = mean(pmet,1); %% Within mean. 
      UDIAG(ligne,4,j) = UDIAG(ligne,4,j) + sum((pmet(:,1)-moyenne).^2)/(n-1);
      UDIAG(ligne,6,j) = UDIAG(ligne,6,j) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
    end
  end
  fprintf('Done! \n');
end
UDIAG(:,[2 4 6],:) = UDIAG(:,[2 4 6],:)/NumberOfBlocks;
disp(' ')
clear pmet temp moyenne CSUP CINF csup cinf n linea iter tmp;    
boxplot = 0;

Info.PlotProperties.CutTop = 0;
Info.PlotProperties.CutBottom = 0;
for j = 1:NumberOfParameters
  [nam,namtex] = get_the_name(j,options_.TeX);
  for i = 1:2:5
    boxplot = boxplot+1;
    eval(['Info.Box' int2str(boxplot) '.Curve1.xdata = xx;'])
    eval(['Info.Box' int2str(boxplot) '.Curve2.xdata = xx;'])
    eval(['Info.Box' int2str(boxplot) '.Curve1.ydata = UDIAG(:,i,j);'])
    eval(['Info.Box' int2str(boxplot) '.Curve2.ydata = UDIAG(:,i+1,j);'])
    eval(['Info.Box' int2str(boxplot) '.Curve1.variablename = ''Pooled moments'';'])
    eval(['Info.Box' int2str(boxplot) '.Curve2.variablename = ''Within moments'';'])
    eval(['Info.Box' int2str(boxplot) '.Curve1.type = ''DiagnosticPooled'' ;'])
    eval(['Info.Box' int2str(boxplot) '.Curve2.type = ''DiagnosticWithin'' ;'])
    if i == 1
      eval(['Info.Box' int2str(boxplot) '.name = '' Interval (' nam ')  '';'])
      if options_.TeX
	eval(['Info.Box' int2str(boxplot) '.texname = '' Interval (' namtex ')  '';'])
      end	
    elseif i == 2
      eval(['Info.Box' int2str(boxplot) '.name = '' m2 (' nam ')  '';'])
      if options_.TeX
	eval(['Info.Box' int2str(boxplot) '.texname = '' \sigma^2 (' namtex ')  '';'])
      end    
    elseif i == 3
      eval(['Info.Box' int2str(boxplot) '.name = '' m3 (' nam ')  '';'])
      if options_.TeX
	eval(['Info.Box' int2str(boxplot) '.texname = '' \sigma^3 (' namtex ')  '';'])
      end
    end
  end
end
% MakeAllFigures()









% $$$ pages = floor(NumberOfParameters/3);
% $$$ k = 0;  
% $$$ for i = 1:pages
% $$$   h = figure('Name','MCMC univariate diagnostic (Brooks and Gelman,1998)');
% $$$   boxplot = 1;
% $$$   if TeX
% $$$     NAMES = [];
% $$$     TEXNAMES = [];
% $$$   end
% $$$   for j = 1:3 % Loop over parameters
% $$$     k = k+1;
% $$$     [nam,namtex] = get_the_name(k,TeX);
% $$$     for crit = 1:3% Loop over criteria
% $$$ 				if crit == 1
% $$$ 	  				plt1 = UDIAG(:,1,k);
% $$$ 	  				plt2 = UDIAG(:,2,k);
% $$$ 	  				namnam  = [nam , ' (Interval)']; 
% $$$ 				elseif crit == 2
% $$$ 	  				plt1 = UDIAG(:,3,k);
% $$$ 	  				plt2 = UDIAG(:,4,k);
% $$$ 	  				namnam  = [nam , ' (m2)'];
% $$$ 				elseif crit == 3    
% $$$ 	  				plt1 = UDIAG(:,5,k);
% $$$ 	  				plt2 = UDIAG(:,6,k);
% $$$ 	  				namnam  = [nam , ' (m3)'];
% $$$ 				end
% $$$ 				if TeX
% $$$ 	  				NAMES = strvcat(NAMES,deblank(namnam));
% $$$ 	  				TEXNAMES = strvcat(TEXNAMES,deblank(namtex));
% $$$ 				end
% $$$ 				subplot(3,3,boxplot);
% $$$ 				plot(xx,plt1,'-b');     % Pooled
% $$$ 				hold on;
% $$$ 				plot(xx,plt2,'-r');     % Within (mean)
% $$$ 				hold off;
% $$$ 				xlim([xx(1) xx(number_of_lines)])
% $$$ 				title(namnam,'Interpreter','none')
% $$$ 				boxplot = boxplot + 1;
% $$$       		end
% $$$     	end
% $$$     	eval(['print -depsc2 ' M_.fname '_udiag' int2str(i)]);
% $$$     	eval(['print -dpdf ' M_.fname '_udiag' int2str(i)]);
% $$$     	saveas(h,[M_.fname '_udiag' int2str(i) '.fig']);
% $$$     	if options_.nograph, close(h), end
% $$$     	if TeX
% $$$       		fprintf(fidTeX,'\\begin{figure}[H]\n');
% $$$       		for jj = 1:size(NAMES,1)
% $$$ 				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
% $$$       		end    
% $$$       		fprintf(fidTeX,'\\centering \n');
% $$$       		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',M_.fname,int2str(i));
% $$$       		fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
% $$$       		fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
% $$$       		fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
% $$$       		fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(i));
% $$$       		fprintf(fidTeX,'\\end{figure}\n');
% $$$       		fprintf(fidTeX,'\n');
% $$$     	end
% $$$ 	end
% $$$   	reste = npar-k;
% $$$   	if reste
% $$$     	if reste == 1
% $$$       		nr = 3;
% $$$       		nc = 1;
% $$$     	elseif reste == 2;
% $$$       		nr = 2;
% $$$       		nc = 3;
% $$$     	end
% $$$     	if TeX
% $$$       		NAMES = [];
% $$$       		TEXNAMES = [];
% $$$     	end
% $$$     	h = figure('Name','MCMC univariate diagnostic (Brooks and Gelman, 1998)');
% $$$     	boxplot = 1;
% $$$     	for j = 1:reste
% $$$       		k = k+1;
% $$$       		[nam,namtex] = get_the_name(k,TeX);
% $$$       		for crit = 1:3
% $$$ 				if crit == 1
% $$$ 	  				plt1 = UDIAG(:,1,k);
% $$$ 	  				plt2 = UDIAG(:,2,k);
% $$$ 	  				namnam  = [nam , ' (Interval)']; 
% $$$ 				elseif crit == 2
% $$$ 	  				plt1 = UDIAG(:,3,k);
% $$$ 	  				plt2 = UDIAG(:,4,k);
% $$$ 	  				namnam  = [nam , ' (m2)'];
% $$$ 				elseif crit == 3    
% $$$ 	  				plt1 = UDIAG(:,5,k);
% $$$ 	  				plt2 = UDIAG(:,6,k);
% $$$ 	  				namnam  = [nam , ' (m3)'];
% $$$ 				end
% $$$ 				if TeX
% $$$ 	  				NAMES = strvcat(NAMES,deblank(namnam));
% $$$ 	  				TEXNAMES = strvcat(TEXNAMES,deblank(namtex));
% $$$ 				end
% $$$ 				subplot(nr,nc,boxplot);
% $$$ 				plot(xx,plt1,'-b');					% Pooled
% $$$ 				hold on;
% $$$ 				plot(xx,plt2,'-r');					% Within (mean)
% $$$ 				hold off;
% $$$ 				xlim([xx(1) xx(number_of_lines)]);
% $$$ 				title(namnam,'Interpreter','none');
% $$$ 				boxplot = boxplot + 1;
% $$$       		end
% $$$     	end
% $$$     	eval(['print -depsc2 ' M_.fname '_udiag' int2str(pages+1)]);
% $$$     	eval(['print -dpdf ' M_.fname '_udiag' int2str(pages+1)]);
% $$$     	saveas(h,[M_.fname '_udiag' int2str(pages+1) '.fig']);
% $$$     	if options_.nograph, close(h), end
% $$$     	if TeX
% $$$       		fprintf(fidTeX,'\\begin{figure}[H]\n');
% $$$       		for jj = 1:size(NAMES,1);
% $$$ 				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
% $$$       		end    
% $$$       		fprintf(fidTeX,'\\centering \n');
% $$$       		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',M_.fname,int2str(pages+1));
% $$$       		if reste == 2
% $$$ 				fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
% $$$ 				fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
% $$$ 				fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
% $$$       		elseif reste == 1
% $$$ 				fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
% $$$ 				fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
% $$$ 				fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
% $$$       		end
% $$$       		fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(pages+1));
% $$$       		fprintf(fidTeX,'\\end{figure}\n');
% $$$       		fprintf(fidTeX,'\n');
% $$$       		fprintf(fidTeX,'% End Of TeX file.');
% $$$       		fclose(fidTeX);
% $$$     	end
% $$$ 	end % if reste > 0
%%
%% Multivariate diagnostic.
%
% $$$   	if TeX
% $$$     	fidTeX = fopen([M_.fname '_MultivariateDiagnostics.TeX'],'w');
% $$$     	fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
% $$$     	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
% $$$     	fprintf(fidTeX,' \n');
% $$$     	NAMES = [];
% $$$   	end

tmp = zeros(NumberOfSimulations*NumberOfBlocks,3);
MDIAG = zeros(NumberOfLines,6);
for b = 1:NumberOfBlocks
  startline = 0;
  for file = 0:LastFileNumber
    eval(['load ' M_.fname '_mh' int2str(file) '_blck' int2str(b)]);
    clear x2 post2;
    tmp((b-1)*NumberOfSimulations+startline+1:(b-1)*NumberOfSimulations+...
	CumulatedNumberOfSimulationsPerFile(file+1),1) = logpo2;
    startline = startline+NumberOfSimulationsPerFile(file+1);
  end	
end
clear logpo2;
tmp(:,2) = kron((1:NumberOfBlocks)',ones(NumberOfSimulations,1));
tmp(:,3) = kron(ones(NumberOfBlocks,1),time'); 
tmp = sortrows(tmp,1);
ligne   = 0;
for iter  = Origin:StepSize:NumberOfSimulations
  ligne = ligne+1;
  linea = ceil(0.5*iter);
  n     = iter-linea+1;
  cinf  = round(n*alpha/2);
  csup  = round(n*(1-alpha/2));
  CINF  = round(NumberOfBlocks*n*alpha/2);
  CSUP  = round(NumberOfBlocks*n*(1-alpha/2));
  temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
  MDIAG(ligne,1) = temp(CSUP,1)-temp(CINF,1);
  moyenne = mean(temp(:,1));%% Pooled mean.
  MDIAG(ligne,3) = sum((temp(:,1)-moyenne).^2)/(NumberOfBlocks*n-1);
  MDIAG(ligne,5) = sum(abs(temp(:,1)-moyenne).^3)/(NumberOfBlocks*n-1);
  for i=1:NumberOfBlocks
    pmet = temp(find(temp(:,2)==i));
    MDIAG(ligne,2) = MDIAG(ligne,2) + pmet(csup,1)-pmet(cinf,1);
    moyenne = mean(pmet,1); %% Within mean. 
    MDIAG(ligne,4) = MDIAG(ligne,4) + sum((pmet(:,1)-moyenne).^2)/(n-1);
    MDIAG(ligne,6) = MDIAG(ligne,6) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
  end
end
MDIAG(:,[2 4 6],:) = MDIAG(:,[2 4 6],:)/NumberOfBlocks;	




FigureProperties.Name = 'MCMC univariate diagnostic (Brooks and Gelman, 1998)';
		    














% $$$  	h = figure('Name','Multivatiate diagnostic');
% $$$   	boxplot = 1;
% $$$   	for crit = 1:3
% $$$     	if crit == 1
% $$$     		plt1 = MDIAG(:,1);
% $$$     		plt2 = MDIAG(:,2);
% $$$     		namnam  = 'Interval'; 
% $$$     	elseif crit == 2
% $$$     		plt1 = MDIAG(:,3);
% $$$     		plt2 = MDIAG(:,4);
% $$$     		namnam  = 'm2';
% $$$     	elseif crit == 3    
% $$$     		plt1 = MDIAG(:,5);
% $$$     		plt2 = MDIAG(:,6);
% $$$     		namnam  = 'm3';
% $$$     	end
% $$$     	if TeX
% $$$     		NAMES = strvcat(NAMES,namnam);
% $$$     	end
% $$$     	subplot(3,1,boxplot);
% $$$     	plot(xx,plt1,'-b');  % Pooled
% $$$     	hold on
% $$$     	plot(xx,plt2,'-r');  % Within (mean)
% $$$     	hold off
% $$$     	xlim([xx(1) xx(number_of_lines)])
% $$$     	title(namnam,'Interpreter','none');
% $$$     	boxplot = boxplot + 1;
% $$$   	end
% $$$   	eval(['print -depsc2 ' M_.fname '_mdiag']);
% $$$   	eval(['print -dpdf ' M_.fname '_mdiag']);
% $$$   	saveas(h,[M_.fname '_mdiag.fig']);
% $$$   	if options_.nograph, close(h), end
% $$$   	if TeX
% $$$     	fprintf(fidTeX,'\\begin{figure}[H]\n');
% $$$     	for jj = 1:3
% $$$       		fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),' ');
% $$$     	end    
% $$$     	fprintf(fidTeX,'\\centering \n');
% $$$     	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_mdiag}\n',M_.fname);
% $$$     	fprintf(fidTeX,'\\caption{Multivariate convergence diagnostics for the Metropolis-Hastings.\n');
% $$$     	fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
% $$$     	fprintf(fidTeX,'the eighty percent interval, the second and third moments. The different \n');
% $$$     	fprintf(fidTeX,'parameters are aggregated using the posterior kernel.}');
% $$$     	fprintf(fidTeX,'\\label{Fig:MultivariateDiagnostics}\n');
% $$$     	fprintf(fidTeX,'\\end{figure}\n');
% $$$     	fprintf(fidTeX,'\n');
% $$$     	fprintf(fidTeX,'% End Of TeX file.');
% $$$     	fclose(fidTeX);
% $$$   	end
% $$$ end 
% $$$ End of if ~options_.nodiagnostic 