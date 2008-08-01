function initvalf(fname)

% function initvalf(fname,varargin)
% reads an initial path from the 'fname' file for exogenous and endogenous variables	
%
% INPUTS
%    fname:         name of the function
%    // period:        period
%    // varargin:      list of arguments following period
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2007 Dynare Team
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

  global M_ oo_ options_
  global y_start_date ex_start_date 

  series = 1;
  if exist(fname) == 2
      eval(fname);
  elseif exist([fname '.xls']) == 2
      [data,names_v]=xlsread([fname '.xls']);
      series = 0;
  elseif exist([fname '.mat']) == 2
      load(fname);
  end
  
% $$$   if length(period) == 2
% $$$     period = dy_date(period(1),period(2));
% $$$   end
% $$$   
% $$$   if period - max(M_.maximum_lag,M_.maximum_lag) < 0
% $$$     error(['INITVALF_: not enough data points in database for number of' ...
% $$$ 	   ' lags. Start later!'])
% $$$   end
% $$$   
% $$$   if nargin > 2
% $$$     if strcmp(upper(varargin{1}),'SERIES')
% $$$       series = 1 ;
% $$$     elseif strcmp(upper(varargin{1}),'MAT')
% $$$       series = 0 ;
% $$$     else
% $$$       error(['INITVALF: unknown option ' varargin{1}])
% $$$     end
% $$$   else
% $$$     series = 0 ;
% $$$   end
% $$$   
% $$$   y1 = floor((period-M_.maximum_lag)/M_.freq);
% $$$   p1 = period-M_.maximum_lag-M_.freq*y1;
% $$$   y_start_date(2) = M_.start_date(2) + p1-1;
% $$$   if y_start_date(2) > M_.freq
% $$$     y_start_date(2) = y_start_date(2) - M_.freq;
% $$$     y1 = y1 + 1;
% $$$   end
% $$$   y_start_date(1) = M_.start_date(1)+y1;
% $$$   
% $$$   y1 = floor((period-M_.maximum_lag)/M_.freq);
% $$$   p1 = period-M_.maximum_lag-M_.freq*y1;
% $$$   ex_start_date(2) = M_.start_date(2) + p1-1;
% $$$   if y_start_date(2) > M_.freq
% $$$     ex_start_date(2) = ex_start_date(2) - M_.freq;
% $$$     y1 = y1 + 1;
% $$$   end
% $$$   ex_start_date(1) = M_.start_date(1)+y1;
% $$$   
% $$$   clear y1, p1;
  
  options_.initval_file = 1;
  oo_.endo_simul = [];
  oo_.exo_simul = [];
  
  for i=1:size(M_.endo_names,1)
    if series == 1
%      x = eval([M_.endo_names(i,:) '(period-M_.maximum_lag:period+options_.periods+M_.maximum_lead-1);']);
      x = eval(M_.endo_names(i,:));
      oo_.endo_simul = [oo_.endo_simul; x'];
    else
      k = strmatch(upper(M_.endo_names(i,:)),names_v,'exact');
      if isempty(k)
	error(['INITVALF: ' M_.endo_names(i,:) ' not found'])
      end
      x = data(:,k);
%      oo_.endo_simul = [oo_.endo_simul; x(period-M_.maximum_lag:period+options_.periods+M_.maximum_lead-1)']; 
      oo_.endo_simul = [oo_.endo_simul; x']; 
    end
  end
  
  for i=1:size(M_.exo_names,1)
    if series == 1
%      x = eval([M_.exo_names(i,:) '(period-M_.maximum_lag:period+options_.periods+M_.maximum_lead-1);']);
      x = eval(M_.exo_names(i,:) );
      oo_.exo_simul = [oo_.exo_simul x];
    else
      k = strmatch(upper(M_.exo_names(i,:)),names_v,'exact');
      if isempty(k)
	error(['INITVALF: ' M_.exo_names(i,:) ' not found'])
      end
      x = data(:,k);
%      oo_.exo_simul = [oo_.exo_simul x(period-M_.maximum_lag:period+options_.periods+M_.maximum_lead-1)]; 
      oo_.exo_simul = [oo_.exo_simul x]; 
    end
  end
    












