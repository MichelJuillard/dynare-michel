% Copyright (C) 2001 Michel Juillard
%
function dynare(fname, varargin)
%	DYNARE ( 'Filename' )
%	This command runs dyanre with specified model file in argument
% 	Filename.
%	The name of model file begins with an alphabetic character, 
%	and has a filename extension of .mod or .dyn.
%	When extension is omitted, a model file with .mod extension
%	is processed.

if ~isstr(fname)
	error ('The argument in DYNARE must be a text string.') ;
end
% Testing if file have extension
% If no extension defalut .mod is added
if isempty(strfind(fname,'.'))
	fname=[fname '.mod'];
% Checking file extension
else
	if ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.MOD') ...
		&& ~strcmp(upper(fname(size(fname,2)-3:size(fname,2))),'.DYN')
		error ('Argument is a file name with .mod or .dyn extension');
	end;
end;
dynareroot = strrep(which('dynare.m'),'dynare.m','');
command = [dynareroot 'dynare_m.exe ' fname] ;
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


% MJ 2/9/99: replace clear function by clear ff_
% MJ 4/7/00: change the path of dynare_m
% MJ 02/26/01: replaced local variable x by fname
% MJ 09/19/01: evaluates mod script in 'base' workspace