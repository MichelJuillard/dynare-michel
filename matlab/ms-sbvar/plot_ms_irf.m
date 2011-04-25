function plot_ms_irf(irf,names,title_)
% function [] = plot_ms_irf(irf,names)
% plots the impulse responses from the output from a ms-sbvar
%
% INPUTS
%   irf should be in the form (percentile x horizon x (nvar x nvar)), if banded otherwise
%     ( horizon x (nvar x nvar) )
%
%   names: character list of the names of the variables
%
%   title: optional super title
%
% The element in position (k,i+j*nvars) of ir is the response of the ith 
% variable to the jth shock at horizon k.  Horizon 0 is the contemporaneous 
% response.

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

  dims = size(irf);
  
  if nargin < 3
    title_ = '';
  end
  
  if (length(dims) == 2)
    % Point IRF (horizon x (nvarsxnvars) )
    horizon = dims(1);
    num_percentiles = 1;
  elseif (length(dims) == 3)
    % Banded IRF
    horizon = dims(2);
    num_percentiles = dims(1);
  else
    error('The impulse response matrix passed to be plotted does not appear to be the correct size');
  end
  
  if size(names,1) ~= nvars
    error('The names passed are not the same length as the number of variables')
  end
  
  if ~exist(fname, 'dir')
      mkdir('.',fname);
  end
  if ~exist([fname '/graphs'])
      mkdir(fname,'graphs');
  end
  
  
  if num_percentiles == 1
    % loop through the shocks
    for s=1:nvars
      shock = zeros(horizon,nvars);
      for i=1:nvars
        shock(:,i) = irf(:,((i-1) + ((s-1)*nvars)+1));
      end
      plot_point_irf_for_shock(shock,nvars,names,names(s,:));
    end
  else
    for s=1:nvars
      shock = zeros(horizon,nvars,num_percentiles);
      for n=1:num_percentiles
        for i=1:nvars
          shock(:,i,n) = irf(n,:,((i-1) + ((s-1)*nvars)+1));
        end
      end
      plot_banded_irf_for_shock(shock,nvars,names,names(s,:));
    end
  end
  
  function [fig] = plot_point_irf_for_shock(irf,nvars,names,shock_name)
    fig = figure('Name',title_);
    for k=1:nvars
      subplot(ceil(sqrt(nvars)), ceil(sqrt(nvars)),k);
      plot(irf(:,k))
      title([names(k,:) ' shock from ' shock_name]);
    end
    %suptitle(title_);
  end
  
  function [fig] = plot_banded_irf_for_shock(irf,nvars, names, shock_name)
    fig = figure('Name',title_);
    npercentiles = size(irf,3);
    for ii=1:nvars
      subplot(ceil(sqrt(nvars)), ceil(sqrt(nvars)),ii);
      for nn=1:npercentiles
        plot(irf(:,ii,nn))
        hold on
      end
      hold off
      title([names(ii,:) ' shock from ' shock_name]);
    end
    %suptitle(title_);
  end
  
end