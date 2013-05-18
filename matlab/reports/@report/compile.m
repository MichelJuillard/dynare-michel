function o = compile(o, varargin)
%function o = compile(o)
% Compile Report Object
%
% INPUTS
%   o            [report]  report object
%   varargin     [char]    allows user to change report compiler for a
%                          given run of compile.
%
% OUTPUTS
%   o     [report]  report object
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

assert(length(varargin) == 0 || length(varargin) == 2, ...
       '@report.compile: calling form: compiler, ''/path/to/compiler''.');
if length(varargin) == 2
    assert(ischar(varargin{1}) && strcmp(lower(varargin{1}), 'compiler'), ...
           '@report.compile: ''compiler'' is the only option.');
    assert(ischar(varargin{2}), ...
           '@report.compile: the argument to ''compiler'' must be a char');
    compiler = varargin{2};
else
    compiler = o.compiler;
end

if ~exist(o.filename, 'file')
    o.write();
end

if isempty(compiler)
    if strncmp(computer, 'MACI', 4) || ~isempty(regexpi(computer, '.*apple.*', 'once'))
        % Add most likely places for pdflatex to exist outside of default $PATH
        [status, compiler] = ...
            system(['PATH=$PATH:/usr/texbin:/usr/local/bin:/usr/local/sbin;' ...
                    'which pdflatex'], '-echo');
    elseif strcmp(computer, 'PCWIN') || strcmp(computer, 'PCWIN64')
        error(['@report.compile: On Windows machines, you must explicitly ' ...
               'provide the ''compiler'' option or set the compiler ' ...
               'variable in the Report class']);
    else % gnu/linux
        [status, compiler] = system('which pdflatex', '-echo');
    end
    assert(status == 0, ...
           '@report.compile: Could not find a tex compiler on your system');
    compiler = strtrim(compiler);
    o.compiler = compiler;
end

if exist('OCTAVE_VERSION')
    status = system([compiler ' ./' o.filename], 0);
else
    status = system([compiler ' ./' o.filename], '-echo');
end
[junk, rfn, junk] = fileparts(o.filename);

if status ~= 0
    error(['@report.compile: There was an error in compiling ' rfn '.pdf.' ...
          '  ' compiler ' returned the error code: ' num2str(status)]);
end
fprintf(1, '\n\nDone.\n')
disp('Your compiled report is located here:');
disp(['     ' pwd filesep rfn '.pdf']);

if ~exist('OCTAVE_VERSION')
    open([pwd filesep rfn '.pdf']);
end
end