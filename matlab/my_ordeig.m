function eval = my_ordeig(t)

% function eval = my_ordeig(t)
% Computes the eigenvalues of a quasi-triangular matrix
%
% INPUTS
%    t:              quasi-triangular matrix
%
% OUTPUTS
%    eval:           eigenvalues
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  n = size(t,2);
  eval = zeros(n,1);
  for i=1:n-1
    if t(i+1,i) == 0
      eval(i) = t(i,i);
    else
      k = i:i+1;
      eval(k) = eig(t(k,k));
      i = i+1;
    end
  end
  if i < n
    eval(n) = t(n,n);
  end
  
      