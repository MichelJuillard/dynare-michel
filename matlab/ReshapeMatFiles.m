function ReshapeMatFiles(type)
  % Rehape and sort (along the mcmc simulations) the mat files generated
  % by DYNARE.
  global M_ options_

  MhDirectoryName = [ CheckPath('metropolis') '/' ];
  
  switch type
   case 'irf'
    IRFfiles = dir([MhDirectoryName M_.fname '_irf*']);
    NumberOfIRFfiles = length(IRFfiles);
    if NumberOfIRFfiles>1
      NumberOfPeriodsPerIRFfiles = ceil(options_.irf/NumberOfIRFfiles);
      reste = options_.irf-NumberOfPeriodsPerIRFfiles*(NumberOfIRFfiles-1);
      idx = 0;
      jdx = 0;
      for f1=1:NumberOfIRFfiles-1
        STOCK_IRF = zeros(NumberOfPeriodsPerIRFfiles,M_.endo_nbr,M_.exo_nbr,B);
        for f2 = 1:NumberOfIRFfiles
          load([MhDirectoryName M_.fname '_irf' int2str(f2)]);
          STOCK_IRF(:,:,:,idx+1:idx+size(stock_irf,4)) = stock_irf(jdx+1:jdx+NumberOfPeriodsPerIRFfiles,:,:,:);
          idx = idx+size(stock_irf,4);
        end
        STOCK_IRF = sort(STOCK_IRF,4);
        save([MhDirectoryName M_.fname '_IRFs' int2str(f1)],'STOCK_IRF');
        jdx = jdx + NumberOfPeriodsPerIRFfiles;
        idx = 0;
      end
      STOCK_IRF = zeros(reste,M_.endo_nbr,M_.exo_nbr,B);
      for f2 = 1:NumberOfIRFfiles
        load([MhDirectoryName M_.fname '_irf' int2str(f2)]);
        STOCK_IRF(:,:,:,idx+1:idx+size(stock_irf,4)) = stock_irf(jdx+1:jdx+reste,:,:,:);
        idx = idx+size(stock_irf,4);
      end
      STOCK_IRF = sort(STOCK_IRF,4);
      save([MhDirectoryName M_.fname '_IRFs' int2str(NumberOfIRFfiles)],'STOCK_IRF');
    else
      load([MhDirectoryName M_.fname '_irf1']);
      STOCK_IRF = sort(stock_irf,4);
      save([MhDirectoryName M_.fname '_IRFs' int2str(1)],'STOCK_IRF');
    end    
    for file = 1:NumberOfIRFfiles
      delete([MhDirectoryName M_.fname '_irf' int2str(file) '.mat'])
    end    
   case { 'smooth' , 'filter' , 'error' , 'innov' , 'forcst' , 'forcst1' }
    TYPEfiles = dir([MhDirectoryName M_.fname '_' type '*']);
    NumberOfTYPEfiles = length(TYPEfiles);
    switch type
     case 'smooth'
      CAPtype  = 'SMOOTH';
      TYPEsize = [ M_.endo_nbr , options_.nobs ];
     case 'filter'
      CAPtype = 'FILTER';
      TYPEsize = [ M_.endo_nbr , options_.nobs + 1 ];% TO BE CHECKED!
     case 'error'
      CAPtype = 'ERROR';
      TYPEsize = [ size(options_.varobs,1) , options_.nobs ];
     case 'innov'
      CAPtype = 'INNOV';
      TYPEsize = [ M_.exo_nbr , options_.nobs ];
     case 'forcst'
      CAPtype = 'forcst';
      TYPEsize = [ M_.endo_nbr , options_.forecast ];
     case 'forcst1'
      CAPtype = 'forcst1';
      TYPEsize = [ M_.endo_nbr , options_.forecast ];
    end
    if NumberOfTYPEfiles>1
      NumberOfPeriodsPerTYPEfiles = ceil( TYPEsize(2)/NumberOfTYPEfiles );
      reste = TYPEsize(2)-NumberOfPeriodsPerTYPEfiles*(NumberOfTYPEfiles-1);
      idx = 0;
      jdx = 0;
      for f1=1:NumberOfIRFfiles-1
        eval(['STOCK_' CAPtype ' = zeros(TYPEsize(1),NumberOfPeriodsPerTYPEfiles,B);'])
        for f2 = 1:NumberOfIRFfiles
          load([MhDirectoryName M_.fname '_' type int2str(f2)]);
          eval(['STOCK_' CAPtype '(:,:,idx+1:idx+size(stock_ ' type ',3))=stock_' type '(:,jdx+1:jdx+NumberOfPeriodsPerTYPEfiles,:);'])
          eval(['idx = idx + size(stock_' type ',3);'])
        end
        eval(['STOCK_' CAPtype ' = sort(STOCK_ ' CAPtype ',3);'])
        save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(f1)],['STOCK_' CAPtype]);
        jdx = jdx + NumberOfPeriodsPerIRFfiles;
        idx = 0;
      end
      eval(['STOCK_' CAPtype ' = zeros(TYPEsize(1),reste,B);'])
      for f2 = 1:NumberOfIRFfiles
        load([MhDirectoryName M_.fname '_' type int2str(f2)]);
        eval(['STOCK_' CAPtype '(:,:,idx+1:idx+size(stock_ ' type ',3))=stock_' type '(:,jdx+1:jdx+reste,:);'])
        eval(['idx = idx + size(stock_' type ',3);'])
      end
      eval(['STOCK_' CAPtype ' = sort(STOCK_ ' CAPtype ',3);'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(NumberOfIRFfiles)],['STOCK_' CAPtype]);
    else
      load([MhDirectoryName M_.fname '_' type '1']);
      eval(['STOCK_' CAPtype ' = sort(stock_' type ',3);'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(1)],['STOCK_' CAPtype ]);      
    end
    for file = 1:NumberOfTYPEfiles
      delete([MhDirectoryName M_.fname '_' type ' int2str(file) '.mat'])
    end
   otherwise
    disp('ReshapeMatFiles :: Unknown format!')
    return
  end