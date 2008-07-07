function eigs = my_ordeig(t)
% function eval = my_ordeig(t)
% Computes the eigenvalues of a quasi-triangular matrix
%
% INPUTS
%    t:              quasi-triangular matrix
%
% OUTPUTS
%    eigs:           eigenvalues
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  n = size(t,2);
  eigs = zeros(n,1);
  i = 1;
  while i <= n
      if i == n
          eigs(n) = t(n,n);
          break;
      elseif t(i+1,i) == 0
          eigs(i) = t(i,i);
      else
          k = i:i+1;
          eigs(k) = eig(t(k,k));
          i = i+1;
      end
      i = i+1;
  end
