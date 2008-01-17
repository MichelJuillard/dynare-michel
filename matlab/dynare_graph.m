function dynare_graph(y,tit,x)

% function dynare_graph(y,tit,x) 
% graphs
%
% INPUT
%   figure_name: name of the figures
%   colors: line colors
%
% OUTPUT
%   none
%
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.

  global dyn_graph

  if nargin < 3
    x = (1:size(y,1))';
  end
  nplot = dyn_graph.plot_nbr + 1; 
  if nplot > dyn_graph.max_nplot
    figure('Name',dyn_graph.figure_name);
    nplot = 1;
  end
  dyn_graph.plot_nbr = nplot;
  subplot(dyn_graph.nr,dyn_graph.nc,nplot);
  
  line_types = dyn_graph.line_types;
  line_type = line_types{1};
  for i=1:size(y,2);
    if length(line_types) > 1
      line_type = line_types{i};
    end
    
    plot(x,y(:,i),line_type);
    hold on
  end
  title(tit);
  hold off