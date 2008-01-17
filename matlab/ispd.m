function test = ispd(A)

% function test = ispd(A)
% Tests if a square matrix is positive definite. 
% 
% INPUTS 
%   o A       [double]   a square matrix. 
%  
% OUTPUTS 
%   o test    [integer]  is equal to one if A is pd, 0 otherwise. 
%
% SPECIAL REQUIREMENTS
%   None.
%  
% part of DYNARE, copyright Dynare Team (2007-2008)
% Gnu Public License.

m = length(A);% I do not test for a square matrix...
test = 1;

for i=1:m
  if ( det( A(1:i, 1:i) ) < 2.0*eps )
    test = 0;
    break
  end
end