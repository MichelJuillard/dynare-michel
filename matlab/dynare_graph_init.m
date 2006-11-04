function dynare_graph_init(figure_name,nplot,line_types,line_width)
% function dynare_graph_init(figure_name,colors) initializes set of
% graphs
% INPUT
%   figure_name: name of the figures
%   colors: line colors
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
  global dyn_graph options_
  
  dyn_graph.fh = figure('Name',figure_name);
  dyn_graph.figure_name = figure_name;
  if nargin > 2
    dyn_graph.line_types = line_types;
  else
    dyn_graph.line_types = options_.graphics.line_types;
  end
  if nargin > 3
    dyn_graph.line_width = line_width;
  else
    dyn_graph.line_width = options_.graphics.line_width;
  end
  
  dyn_graph.plot_nbr = 0;
  
  switch(nplot)
   case 1
    dyn_graph.nc = 1;
    dyn_graph.nr = 1;
   case 2
    dyn_graph.nc = 1;
    dyn_graph.nr = 2;
   case 3
    dyn_graph.nc = 1;
    dyn_graph.nr = 3;
   case 4
    dyn_graph.nc = 2;
    dyn_graph.nr = 2;
   case 5
    dyn_graph.nc = 3;
    dyn_graph.nr = 2;
   case 6
    dyn_graph.nc = 3;
    dyn_graph.nr = 2;
   case 7
    dyn_graph.nc = 4;
    dyn_graph.nr = 2;
   case 8
    dyn_graph.nc = 4;
    dyn_graph.nr = 2;
   otherwise
    dyn_graph.nc = min(nplot,options_.graphics.ncols);
    dyn_graph.nr = min(ceil(nplot/dyn_graph.nc),options_.graphics.nrows);
  end
  dyn_graph.max_nplot = dyn_graph.nc*dyn_graph.nr;
    