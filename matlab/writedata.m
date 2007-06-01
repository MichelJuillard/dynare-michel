function writedata(fname)
% function writedata(fname)
% store endogenous and exogenous variables in a XLS spreadsheet file 
% INPUT
%   fname: name of the XLS file
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright (2007)
% Gnu Public License.
  global M_ oo_
  S=[fname '_endo.xls'];
  n=size(oo_.endo_simul,2);
  delete(S);
  S=upper(cellstr(M_.endo_names));
  S1=cellstr([num2str((1:n)')  char(65*ones(1,n))']);
  xlswrite([fname '_endo'], S', 'endogenous', 'B1');
  xlswrite([fname '_endo'], S1, 'endogenous', 'A2');
  xlswrite([fname '_endo'], oo_.endo_simul', 'endogenous', 'B2');
  S=[fname '_exo.xls'];
  n=size(oo_.exo_simul,1);
  delete(S);
  S=upper(cellstr(M_.exo_names));
  S1=cellstr([num2str((1:n)')  char(65*ones(1,n))']);
  xlswrite([fname '_exo'], S','exogenous', 'B1');
  xlswrite([fname '_exo'], S1, 'exogenous', 'A2');
  xlswrite([fname '_exo'], oo_.exo_simul,'exogenous', 'B2');
return;