function []=dynare(fname)
// Copyright (C) 2001 Michel Juillard
// 
// DYNARE ( 'Filename' )
// This command runs the .MOD file specified in Filename.
// Filename could be enter with or without the .MOD extension.
 
if ~(type(fname)==10) then
  error('The argument in DYNARE must be a text string.');
end
k = strindex(fname,'.')
if k then
  fname = part(fname,1:k-1);
end

command = 'e:/dynare/scilab/dynare_s '+fname;
 
disp(command);

unix_w(command);

exec(fname+'.sci',-1) ;
 
endfunction 
// MJ 2/9/99: replace clear function by clear ff_
// MJ 4/7/00: change the path of dynare_m
// MJ 02/26/01: replaced local variable x by fname
// MJ 09/14/01: Scilab translation
