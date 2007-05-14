function ReshapeMatFiles(type, type2)
  % Reshape and sort (along the mcmc simulations) the mat files generated
  % by DYNARE.
  %
  % 4D-arrays are splitted along the first dimension.
  % 3D-arrays are splitted along the second dimension.
  %
  %   
  global M_ options_

if nargin==1, 
  MhDirectoryName = [ CheckPath('metropolis') '/' ];
else
  if strcmpi(type2,'posterior')
  MhDirectoryName = [CheckPath('metropolis') '/' ];
elseif strcmpi(type2,'gsa')
if options_.opt_gsa.morris==1,
  MhDirectoryName = [CheckPath('GSA\SCREEN') '/' ];
elseif options_.opt_gsa.morris==2,
  MhDirectoryName = [CheckPath('GSA\IDENTIF') '/' ];
else
  MhDirectoryName = [CheckPath('GSA') '/' ];
end
else
  MhDirectoryName = [CheckPath('prior') '/' ];
end  
end
  switch type
   case 'irf_dsge'
    CAPtype  = 'IRF_DSGE';
    TYPEsize = [ options_.irf , M_.endo_nbr , M_.exo_nbr ];
    TYPEarray = 4;    
   case 'irf_bvardsge'
    CAPtype  = 'IRF_BVARDSGE';
    TYPEsize = [ options_.irf , size(options_.varobs,1) , M_.exo_nbr ];
    TYPEarray = 4;      
   case 'smooth'
    CAPtype  = 'SMOOTH';
    TYPEsize = [ M_.endo_nbr , options_.nobs ];
    TYPEarray = 3;
   case 'filter'
    CAPtype = 'FILTER';
    TYPEsize = [ M_.endo_nbr , options_.nobs + 1 ];% TO BE CHECKED!
    TYPEarray = 3;
   case 'error'
    CAPtype = 'ERROR';
    TYPEsize = [ size(options_.varobs,1) , options_.nobs ];
    TYPEarray = 3;
   case 'innov'
    CAPtype = 'INNOV';
    TYPEsize = [ M_.exo_nbr , options_.nobs ];
    TYPEarray = 3;
   case 'forcst'
    CAPtype = 'FORCST';
    TYPEsize = [ M_.endo_nbr , options_.forecast ];
    TYPEarray = 3;
   case 'forcst1'
    CAPtype = 'FORCST1';
    TYPEsize = [ M_.endo_nbr , options_.forecast ];
    TYPEarray = 3;
   otherwise
    disp('ReshapeMatFiles :: Unknown argument!')
    return
  end
  
  TYPEfiles = dir([MhDirectoryName M_.fname '_' type '*.mat']);
  NumberOfTYPEfiles = length(TYPEfiles);
  B = options_.B;
  
  switch TYPEarray
   case 4
    if NumberOfTYPEfiles > 1
      NumberOfPeriodsPerTYPEfiles = ceil(TYPEsize(1)/NumberOfTYPEfiles);
      foffset = NumberOfTYPEfiles-floor(TYPEsize(1)/NumberOfPeriodsPerTYPEfiles);
      reste = TYPEsize(1)-NumberOfPeriodsPerTYPEfiles*(NumberOfTYPEfiles-foffset);
      idx = 0;
      jdx = 0;
      for f1=1:NumberOfTYPEfiles-foffset
        eval(['STOCK_' CAPtype ' = zeros(NumberOfPeriodsPerTYPEfiles,TYPEsize(2),TYPEsize(3),B);'])
        for f2 = 1:NumberOfTYPEfiles
          load([MhDirectoryName M_.fname '_' type int2str(f2)]);
          eval(['STOCK_' CAPtype '(:,:,:,idx+1:idx+size(stock_' type ',4))=stock_' type '(jdx+1:jdx+NumberOfPeriodsPerTYPEfiles,:,:,:);'])
          eval(['idx = idx + size(stock_' type ',4);'])          
        end
        %eval(['STOCK_' CAPtype ' = sort(STOCK_' CAPtype ',4);'])
        save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(f1)],['STOCK_' CAPtype]);
        jdx = jdx + NumberOfPeriodsPerTYPEfiles;
        idx = 0;
      end
      if reste
      eval(['STOCK_' CAPtype ' = zeros(reste,TYPEsize(2),TYPEsize(3),B);'])
      for f2 = 1:NumberOfTYPEfiles
        load([MhDirectoryName M_.fname '_' type int2str(f2)]);
        eval(['STOCK_' CAPtype '(:,:,:,idx+1:idx+size(stock_' type ',4))=stock_' type '(jdx+1:jdx+reste,:,:,:);'])
        eval(['idx = idx + size(stock_' type ',4);'])
      end
      %eval(['STOCK_' CAPtype ' = sort(STOCK_' CAPtype ',4);'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(NumberOfTYPEfiles-foffset+1)],['STOCK_' CAPtype]);  
      end
    else
      load([MhDirectoryName M_.fname '_' type '1']);
      %eval(['STOCK_' CAPtype ' = sort(stock_' type ',4);'])
      eval(['STOCK_' CAPtype ' = stock_' type ';'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(1)],['STOCK_' CAPtype ]);
    end
    % Original file format may be useful in some cases...
    % for file = 1:NumberOfTYPEfiles
    %  delete([MhDirectoryName M_.fname '_' type int2str(file) '.mat'])
    % end
   case 3
    if NumberOfTYPEfiles>1
      NumberOfPeriodsPerTYPEfiles = ceil( TYPEsize(2)/NumberOfTYPEfiles );
      reste = TYPEsize(2)-NumberOfPeriodsPerTYPEfiles*(NumberOfTYPEfiles-1);
      idx = 0;
      jdx = 0;
      for f1=1:NumberOfTYPEfiles-1
        eval(['STOCK_' CAPtype ' = zeros(TYPEsize(1),NumberOfPeriodsPerTYPEfiles,B);'])
        for f2 = 1:NumberOfTYPEfiles
          load([MhDirectoryName M_.fname '_' type int2str(f2)]);
          eval(['STOCK_' CAPtype '(:,:,idx+1:idx+size(stock_ ' type ',3))=stock_' type '(:,jdx+1:jdx+NumberOfPeriodsPerTYPEfiles,:);'])
          eval(['idx = idx + size(stock_' type ',3);'])
        end
        %eval(['STOCK_' CAPtype ' = sort(STOCK_' CAPtype ',3);'])
        save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(f1)],['STOCK_' CAPtype]);
        jdx = jdx + NumberOfPeriodsPerTYPEfiles;
        idx = 0;
      end
      eval(['STOCK_' CAPtype ' = zeros(TYPEsize(1),reste,B);'])
      for f2 = 1:NumberOfTYPEfiles
        load([MhDirectoryName M_.fname '_' type int2str(f2)]);
        eval(['STOCK_' CAPtype '(:,:,idx+1:idx+size(stock_' type ',3))=stock_' type '(:,jdx+1:jdx+reste,:);'])
        eval(['idx = idx + size(stock_' type ',3);'])
      end
      %eval(['STOCK_' CAPtype ' = sort(STOCK_' CAPtype ',3);'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(NumberOfTYPEfiles)],['STOCK_' CAPtype]);
    else
      load([MhDirectoryName M_.fname '_' type '1']);
      %eval(['STOCK_' CAPtype ' = sort(stock_' type ',3);'])
      eval(['STOCK_' CAPtype ' = stock_' type ';'])
      save([MhDirectoryName M_.fname '_' CAPtype 's' int2str(1)],['STOCK_' CAPtype ]);      
    end
    % Original file format may be useful in some cases...
    % for file = 1:NumberOfTYPEfiles
    %   delete([MhDirectoryName M_.fname '_' type  int2str(file) '.mat'])
    % end
  end