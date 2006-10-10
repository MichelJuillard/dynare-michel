function SampleAddress = selec_posterior_draws(SampleSize,info,filepath,filename)
% Selects a sample of draws from the posterior distribution.
% 
% INPUTS
%   o SampleSize     [integer]  Size of the sample to build
%   o info           [integer]  Ff 1 then posterior draws are saved on disk. 
% OUTPUTS 
%   o SampleAddress  [integer]  A (SampleSize*4) array, each line specifies the 
%                               location of a posterior draw: 
%                                  Column 2 --> Chain number
%                                  Column 3 --> (mh) File number    
%                                  Column 4 --> (mh) line number
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
% 
% SPECIAL THANKS 
%   o Antoine.   
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
  global M_ options_ estim_params_
    
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
  
  SampleAddress = zeros(SampleSize,4);
  for i = 1:SampleSize
    ChainNumber = ceil(rand*mh_nblck);
    DrawNumber  = ceil(rand*NumberOfDraws);
    SampleAddress(i,1) = DrawNumber;
    SampleAddress(i,2) = ChainNumber;
    if DrawNumber <= MAX_nruns-FirstLine+1
      MhFileNumber = FirstMhFile;
      MhLineNumber = FirstLine+DrawNumber-1;
    else
      DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
      MhFileNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns); 
      MhLineNumber = DrawNumber-(MhFileNumber-FirstMhFile-1)*MAX_nruns;
    end
    SampleAddress(i,3) = MhFileNumber;
    SampleAddress(i,4) = MhLineNumber;
  end
  SampleAddress = sortrows(SampleAddress,[3 2]);
  if nargin == 2 & info
    if SampleSize <= MAX_nruns% The posterior draws are saved in one file.
      pdraws = zeros(SampleSize,npar);
      old_mhfile = 0;
      old_mhblck = 0;
      for i = 1:SampleSize
        mhfile = SampleAddress(i,3);
        mhblck = SampleAddress(i,2);
        if (mhfile ~= old_mhfile) | (mhblck ~= old_mhblck)
          load([fname '_mh' num2str(mhfile) '_blck' num2str(mhblck) '.mat'],'x2')
        end
        pdraws(i,:) = x2(SampleAddress(i,4),:);
        old_mhfile = mhfile;
        old_mhblck = mhblck;
      end
      clear('x2')
      save([fname '_posterior_draws'],'pdraws')
    else% The posterior draws are saved in xx files.
      NumberOfFiles = ceil(SampleSize/MAX_nruns);
      NumberOfLines = SampleSize - (NumberOfFiles-1)*MAX_nruns;
      linee = 0;
      fnum  = 1;
      pdraws = zeros(MAX_nruns,npar);
      old_mhfile = 0;
      old_mhblck = 0;
      for i=1:SampleSize
        linee = linee+1;
        mhfile = SampleAddress(i,3);
        mhblck = SampleAddress(i,2);
        if (mhfile ~= old_mhfile) | (mhblck ~= old_mhblck)
          load([fname '_mh' num2str(mhfile) '_blck' num2str(mhblck) '.mat'],'x2')
        end
        pdraws(linee,:) = x2(SampleAddress(i,4),:);
        old_mhfile = mhfile;
        old_mhblck = mhblck;
        if fnum < NumberOfFiles && linee == MAX_nruns
          linee = 0;
          save([fname '_posterior_draws' num2str(fnum)],'pdraws')
          fnum = fnum+1;
          if fnum < NumberOfFiles
            pdraws = zeros(MAX_nruns,npar);
          else
            pdraws = zeros(NumberOfLines,npar);
          end
        end
      end
      save([fname '_posterior_draws' num2str(fnum)],'pdraws')
    end    
  end