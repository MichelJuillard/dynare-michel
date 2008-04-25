function E = gensylv(fake,A,B,C,D)
% Solves a Sylvester equation.
%
% INPUTS
%   fake     Unused argument (for compatibility with the mex file)
%   A
%   B
%   C
%   D
%    
% OUTPUTS
%   E
%    
% ALGORITHM
%   none.
%
% SPECIAL REQUIREMENTS
%   none.
%  
%  
% part of DYNARE, copyright Dynare Team (1996-2008)
% Gnu Public License.
C  = kron(C,C); 
x0 = sylvester3(A,B,C,D);
E  = sylvester3a(x0,A,B,C,D);