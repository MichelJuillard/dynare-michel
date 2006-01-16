function Draws = GetAllPosteriorDraws(column,FirstMhFile,FirstLine,TotalNumberOfMhFile,NumberOfDraws)
% stephane.adjemian@ens.fr [09-09-2005]
global M_ options_

nblck = options_.mh_nblck; 
iline = FirstLine; 
linee = 1;
DirectoryName = CheckPath('metropolis');
Draws = zeros(NumberOfDraws*nblck,1);


for file = FirstMhFile:TotalNumberOfMhFile
  for blck = 1:nblck
    load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'x2')
    NumberOfLines = size(x2(iline:end,:),1);
    Draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
    linee = linee+NumberOfLines;
  end
  iline = 1;
end