function r = rcond(A)
% Computes reciprocal condition number
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
  r = 1/(norm(A,1) * norm(inv(A), 1));
