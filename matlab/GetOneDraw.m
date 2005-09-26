function xparams = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,nn,DirectoryName)
% stephane.adjemian@ens.fr [09-25-2005]
global options_ M_

ChainNumber = ceil(rand*options_.mh_nblck);
DrawNumber  = ceil(rand*NumberOfDraws);

if DrawNumber <= nn-FirstLine+1
  MhFilNumber = FirstMhFile;
  MhLine = FirstLine+DrawNumber-1;
else
  DrawNumber  = DrawNumber-(nn-FirstLine+1);
  MhFilNumber = FirstMhFile+ceil(DrawNumber/nn); 
  MhLine = DrawNumber-(MhFilNumber-FirstMhFile-1)*nn;
end

load([DirectoryName '\' M_.fname '_mh' int2str(MhFilNumber) '_blck' int2str(ChainNumber) '.mat' ],'x2');
xparams = x2(MhLine,:);