function o = createGraph(o)
%function o = createGraph(o)
% Create the graph
%
% INPUTS
%   o   - Graph Object
%
% OUTPUTS
%   o   - Graph Object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013 Dynare Team
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

assert(~isempty(o.data));
assert(isa(o.data, 'dynSeries')) ;

if ~isempty(o.figname)
    warning('Will overwrite %s with new graph\n', o.figname);
end

%o = readConfig(o);

disp('creating plot..........');
h = figure('visible','off');
hold on;
box on;
if o.grid
    grid on;
    set(gca, 'GridLineStyle', '--');
end
%set(0, 'CurrentFigure',h);
%set(h, 'PaperPositionMode', 'auto');
%set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);

if strcmpi(o.seriestoplot, 'all')
    data = o.data.data;
else
    data = o.data{o.seriestoplot{:}}.data;
end

x=[1:1:o.data.nobs];
xlabels=getDatesCellStringArray(o.data.time);

plot(x, data);

if ~isempty(o.shade)
    x1 = strmatch(lower(o.shade{1}), xlabels, 'exact');
    x2 = strmatch(lower(o.shade{2}), xlabels, 'exact');
    yrange = get(gca, 'YLim');

    % From ShadePlotForEmpahsis (Matlab Exchange)
    % use patch bc area doesn't work with matlab2tikz
    patch([repmat(x1, 1, 2) repmat(x2, 1, 2)], [yrange fliplr(yrange)], ...
          'b', 'FaceAlpha', .2);
end

set(gca,'XTick', x);
set(gca,'XTickLabel', xlabels);

if o.legend
    if strcmpi(o.seriestoplot, 'all')
        lh = legend(o.data.name);
    else
        lh = legend(o.seriestoplot{:});
    end
    set(lh, 'orientation', o.legend_orientation);
    set(lh, 'Location', o.legend_location);
    set(lh, 'FontSize', o.legend_font_size);
    legend('boxoff');
end

if ~isempty(o.xlabel)
    xlabel(['$\textbf{\footnotesize ' o.xlabel '}$'],'Interpreter','LaTex');
end

if ~isempty(o.ylabel)
    ylabel(['$\textbf{\footnotesize ' o.ylabel '}$'],'Interpreter','LaTex');
end

if ~isempty(o.title)
    title( o.title , 'interpreter', 'none', 'FontSize', 20);
end
drawnow;

o.figname = ['figure-' num2str(cputime) '.tex'];
disp('  converting to tex....');
matlab2tikz('filename', o.figname, ...
            'showInfo', false, ...
            'showWarnings', false, ...
            'checkForUpdates', false);

grid off;
box off;
hold off;
close(h);
clear h;
end
