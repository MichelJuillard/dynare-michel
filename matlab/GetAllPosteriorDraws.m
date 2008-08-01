function Draws = GetAllPosteriorDraws(column,FirstMhFile,FirstLine,TotalNumberOfMhFile,NumberOfDraws)

% function Draws = GetAllPosteriorDraws(column,FirstMhFile,FirstLine,TotalNumberOfMhFile,NumberOfDraws)
% Gets all posterior draws
%
% INPUTS
%    column:               column
%    FirstMhFile:          first mh file 
%    FirstLine:            first line
%    TotalNumberOfMhFile:  total number of mh file 
%    NumberOfDraws:        number of draws

% OUTPUTS
%    Draws:                draws from posterior distribution
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2008 Dynare Team
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