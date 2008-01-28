function read_data_

% function read_data_
% reads endogenous and exogenous variables from a text file 
% Used by datafile option in simulate
%
% INPUT
%   none
%
% OUTPUT
%   none
%
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright Dynare Team (2007-2008)
% Gnu Public License.

  global options_ M_ oo_;
  dname= options_.datafile;
  fid = fopen([dname '_endo.dat'],'r');
  names_line = fgetl(fid);
  allVariables = '';
  positions = ones(0);
  while (any(names_line))
    [chopped,names_line] = strtok(names_line);
    allVariables = strvcat(allVariables, chopped);
    positions = [positions ; strmatch(chopped,M_.endo_names,'exact')];
  end
  Values=fscanf(fid,'%f',inf);
  Values=reshape(Values,M_.endo_nbr,size(Values,1)/M_.endo_nbr);
  oo_.endo_simul=Values(positions,:);
  fclose(fid);
  
  fid = fopen([dname '_exo.dat'],'r');
  names_line = fgetl(fid);
  allVariables = '';
  positions = ones(0);
  while (any(names_line))
    [chopped,names_line] = strtok(names_line);
    allVariables = strvcat(allVariables, chopped);
    positions = [positions ; strmatch(chopped,M_.exo_names,'exact')];
  end
  Values=fscanf(fid,'%f',inf);
  Values=reshape(Values,M_.exo_nbr,size(Values,1)/M_.exo_nbr);
  oo_.exo_simul=(Values(positions,:))';
  fclose(fid);
  %disp([allVariables M_.endo_names]);
  %disp(positions);
  
end