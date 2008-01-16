
function dynare(fname, varargin)

% function dynare(fname, varargin)
%	This command runs dynare with specified model file in argument
% 	Filename.
%	The name of model file begins with an alphabetic character, 
%	and has a filename extension of .mod or .dyn.
%	When extension is omitted, a model file with .mod extension
%	is processed.
%
% INPUTS
%   fname:      file name
%   varargin:   list of arguments following fname
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none
%  
% part of DYNARE, copyright Dynare Team (2001-2008)
% Gnu Public License.


MATLAB  = ver('matlab');
% FIXME:
% It's not satisfactory to convert string versions into numbers, and to
% compare these numbers:
% - conversion will fail if version = 1.2.3
% - it will give 7.10 < 7.9
VERSION = str2num(MATLAB.Version);

dynareroot = strrep(which('dynare.m'),'dynare.m','');

if (VERSION <= 7.4)
    addpath([dynareroot 'mex/2007a/'])
elseif (VERSION == 7.5)
    addpath([dynareroot 'mex/2007b/'])
end

if ~isstr(fname)
  error ('The argument in DYNARE must be a text string.') ;
end
% Testing if file have extension
% If no extension defalut .mod is added
if isempty(strfind(fname,'.'))
  fname1 = [fname '.dyn'];
  d = dir(fname1);
  if length(d) == 0
    fname1 = [fname '.mod'];
  end
  fname = fname1;
  % Checking file extension
else
  if ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.MOD') ...
	&& ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.DYN')
    error ('Argument is a file name with .mod or .dyn extension');
  end;
end;
d = dir(fname);
if length(d) == 0
  disp(['DYNARE: can''t open ' fname])
  return
end

command = ['"' dynareroot 'dynare_m" ' fname] ;
for i=2:nargin
  command = [command ' ' varargin{i-1}];
end
[status, result] = system(command);
if status
  error(result)
end

if ~ isempty(find(abs(fname) == 46))
	fname = fname(:,1:find(abs(fname) == 46)-1) ;
end
evalin('base',fname) ;