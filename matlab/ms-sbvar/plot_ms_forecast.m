function plot_ms_forecast(forecast,title_)
% function [] = plot_ms_forecast(forecast,names)
% plots the forecast from the output from a ms-sbvar
%
% INPUTS
%   forecast should be in the form (percentile x horizon x nvar ), if banded otherwise
%     ( horizon x nvar )
%
%   title: optional super title
%

% Copyright (C) 2011 Dynare Team
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

  nc = 2;
  nr = 2;
  nvars = M_.endo_nbr;
  endo_names = M_.endo_names;
  fname = M_.fname;
  
  var_list = endo_names(1:M_.orig_endo_nbr,:);

  i_var = [];
  for i = 1:size(var_list)
      tmp = strmatch(var_list(i,:),endo_names,'exact');
      if isempty(tmp)
          error([var_list(i,:) ' isn''t and endogenous variable'])
      end
      i_var = [i_var; tmp];
  end
  nvar = length(i_var);
  
  dims = size(forecast);
  
  if nargin < 2
    title_ = '';
  end
  
  if (length(dims) == 2)
    % Point Forecast (horizon x nvars )
    horizon = dims(1);
    num_percentiles = 1;
  elseif (length(dims) == 3)
    % Banded Forecast
    horizon = dims(2);
    num_percentiles = dims(1);  
  else
    error('The impulse response matrix passed to be plotted does not appear to be the correct size');
  end

  if ~exist(fname, 'dir')
      mkdir('.',fname);
  end
  if ~exist([fname '/graphs'])
      mkdir(fname,'graphs');
  end

  n_fig = 1;
  if num_percentiles == 1
    plot_point_forecast(forecast,nvars,endo_names);
    
  else
    plot_banded_forecast(forecast,nvars,endo_names,num_percentiles);
  end

  function plot_point_forecast(forecast,nvars,names)
    fig = figure('Name','Forecast (I)'); 
    m = 1;
    for j=1:nvars
      if m > nr*nc
        %eval(['print -deps ' fname '/graphs/forcst' int2str(n_fig) '.eps;'])
        n_fig =n_fig+1;
        eval(['figure(''Name'',''Forecast (' int2str(n_fig) ')'');']);
        m = 1;
      end
      subplot(nr,nc,m);
      vn = deblank(names(i_var(j),:));
      plot(forecast(:,j))
      title(vn,'Interpreter','none');
      %suptitle(title_);
      grid on;
      m = m+1;
    end
    if m > 1
        %eval(['print -deps ' fname '/graphs/forcst' int2str(n_fig) '.eps;'])
    end
  end
  
  function plot_banded_forecast(forecast,nvars,names,num_percentiles)
    fig = figure('Name',[title_ ' (1)']); 
    %suptitle(title_);
    m = 1;
    
    for j=1:nvars
      if m > nr*nc
        suptitle(title_);
        %eval(['print -deps ' fname '/graphs/forcst' int2str(n_fig) '.eps;'])
        n_fig =n_fig+1;
        eval(['figure(''Name'',''' title_ ' (' int2str(n_fig) ')'');']);
        m = 1;
      end
      subplot(nr,nc,m);
      vn = deblank(names(i_var(j),:));
      for k=1:num_percentiles
        if ceil(num_percentiles/2) == k
          plot(forecast(k,:,j),'LineWidth',1.5)
        else
          plot(forecast(k,:,j),'LineWidth',1.1)
        end
        if k==1
          hold on;
        end
      end
      title(vn,'Interpreter','none');
      hold off
      grid on;
      m = m+1;
    end
    if m > 1
        %eval(['print -deps ' fname '/graphs/forcst' int2str(n_fig) '.eps;'])
    end
  end
  

end