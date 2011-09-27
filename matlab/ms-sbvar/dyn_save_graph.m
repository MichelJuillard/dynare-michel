function dyn_save_graph(dirname,graph_name,graph_formats,TeX,names,texnames,caption)
% function dyn_graph_save(graph_name,graph_formats,TeX)
% saves Dynare graphs
%
% INPUTS
%   graph_name    (string)  name of the graph (used as file name)
%   graph_formats (struct)  list of graph formats to be used
%   TeX           (logical) whether to make TeX snippet
%
% OUTPUTS
%   none
%    
% SPECIAL REQUIREMENTS
%    none
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

    create_dir(dirname);
    graph_name = [dirname filesep regexprep(graph_name,' ','_')];
    if nargin <= 2
        TeX = 0;
    elseif nargin <= 4
        names = {};
        texnames = {};
    elseif nargin <= 6
        caption = '';
    end
    
    if graph_formats.eps || TeX
        print([ graph_name '.eps' ],'-depsc2');
    end
    if graph_formats.pdf && ~exist(OCTAVE_VERSION)
        print(graph_name,'-dpdf');
    end
    if graph_formats.fig && ~exist(OCTAVE_VERSION)
        print(graph_name,'-dfig');
    end
    
    if TeX
        fh = fopen([graph_name '.tex'],'w');
        for i=1:length(names)
            fprintf(fh,'\\psfrag{%s}[1][][0.5][0]{%s}\n',names{i},texnames{i});
        end
        fprintf(fh,'\\centering \n');
        fprintf(fh,'\\includegraphics[scale=0.5]{%s}\n',graph_name);
        if caption
            fprintf(fh,'\\caption{%s}',caption);
        end
        fprintf(fh,'\\label{Fig:%s}\n',graph_name);
        fprintf(fh,'\\end{figure}\n');
        fprintf(fh,'\n');
        fprintf(fh,'%% End of TeX file.\n');
        fclose(fh);
    end
