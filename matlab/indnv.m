function a=indnv(x,y)

% function a=indnv(x,y)
% Finds the elements indices of one vector in an other one
%
% INPUTS
%    x:         column vector
%    y:         column vector
%
% OUTPUTS
%    a:         vector of elements position of x in y
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.

a = zeros(size(x)) ;

for i = 1:size(x,1)
  j = find( x(i) == y );
  if isempty(j)
    a(i) = NaN;
  else
    a(i) = j;
  end
end



