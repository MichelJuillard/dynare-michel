function Draws = GetAllPosteriorDraws(column,FirstMhFile,FirstLine,TotalNumberOfMhFile,NumberOfDraws)
% stephane.adjemian@ens.fr [09-09-2005]
global M_ options_

nblck = options_.mh_nblck; 
iline = FirstLine; 
linee = 1;
DirectoryName = CheckPath('metropolis');
Draws = zeros(NumberOfDraws*nblck,1);
logpo = zeros(NumberOfDraws*nblck,1);
ipost=0;
if column<=0, 
  column=1;  
  ipost=1;
end
iline0=iline;
for blck = 1:nblck
  iline=iline0;
  for file = FirstMhFile:TotalNumberOfMhFile
    load([DirectoryName '/'  M_.fname '_mh' int2str(file) '_blck' int2str(blck)],'x2','logpo2')
    NumberOfLines = size(x2(iline:end,:),1);
    Draws(linee:linee+NumberOfLines-1) = x2(iline:end,column);
    logpo(linee:linee+NumberOfLines-1) = logpo2(iline:end);
    linee = linee+NumberOfLines;
    iline = 1;
  end
end

if ipost,
  Draws=logpo;
end