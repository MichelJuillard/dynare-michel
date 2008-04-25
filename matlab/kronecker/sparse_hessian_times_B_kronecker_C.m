function D = sparse_hessian_times_B_kronecker_C(A,B,C)
% Computes A * kron(B,C) where A is a sparse matrix.
%
% INPUTS
%   A  [double] mA*nA matrix.
%   B  [double] mB*nB matrix.
%   C  [double] mC*nC matrix.
%  
% OUTPUTS
%   D  [double] mA*(nC*nB) or mA*(nB*nB) matrix.
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

D = A_times_B_kronecker_C(A,B,C);