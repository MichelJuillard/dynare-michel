function [xparams, logpost]=metropolis_draw(init)

% function [xparams, logpost]=metropolis_draw(init) 
% Builds draws from metropolis
%
% INPUTS:
%   init:              scalar equal to 1 (first call) or 0
%
% OUTPUTS:
%   xparams:           vector of estimated parameters
%   logpost:           log of posterior density
%   
% SPECIAL REQUIREMENTS
%   none
%  
% part of DYNARE, copyright Dynare Team (2003-2007)
% Gnu Public License.
  



  global options_ estim_params_ M_
  persistent mh_nblck NumberOfDraws fname FirstLine FirstMhFile MAX_nruns
  
  if init
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    MhDirectoryName = CheckPath('metropolis');
    fname = [ MhDirectoryName '/' M_.fname];
    load([ fname '_mh_history']);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine; 
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); 
    LastMhFile = TotalNumberOfMhFiles; 
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
    mh_nblck = options_.mh_nblck;
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

  load( [ fname '_mh' int2str(MhFilNumber) '_blck' int2str(ChainNumber) '.mat' ],'x2','logpo2');
  xparams = x2(MhLine,:);
  logpost= logpo2(MhLine);