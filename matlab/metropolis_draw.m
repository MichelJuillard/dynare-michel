function xparam=metropolis_draw(init)
  global options_ estim_params_
  persistent mh_nblck NumberofDraws fname FirstLine FirstMhFile MAX_nruns
  
  if init
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    DirectoryName = CheckPath('metropolis');
    fname = [ DirectoryName '/' M_.fname]
    load([ fname '_mh_history'])
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine; 
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); LastMhFile = TotalNumberOfMhFiles; 
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop* ...
					       TotalNumberOfMhDraws);
    MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
    mh_nblck = options_.nblck;
    return
  end
  
  ChainNumber = ceil(rand*mh_nblck);
  DrawNumber  = ceil(rand*NumberOfDraws);

  if DrawNumber <= MAX_nruns-FirstLine+1
    MhFilNumber = FirstMhFile;
    MhLine = FirstLine+DrawNumber-1;
  else
    DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
    MhFilNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns); 
    MhLine = DrawNumber-(MhFilNumber-FirstMhFile-1)*MAX_nruns;
  end

  load( [ fname '_mh' int2str(MhFilNumber) '_blck' int2str(ChainNumber) '.mat' ],'x2');
  xparams = x2(MhLine,:);
