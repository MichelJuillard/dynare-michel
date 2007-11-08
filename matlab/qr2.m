function [Q,R] = qr2(X)
% This routine performs a qr decomposition of matrix X such that the
% diagonal scalars of the upper-triangular matrix R are positive. If X 
% is a full (column) rank matrix, then R is also the cholesky
% factorization of X'X. This property is needed for the Del Negro 
% & Schorfheides's identification scheme.
% 
% INPUTS 
%   See matlab's documentation.
%     
% OUTPUTS 
%   See matlab's documentation.
%
% ALGORITHM
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% stephane.adjemian@ens.fr [12-07-2005]
% part of DYNARE, copyright Dynare Team (2007)
% Gnu Public License.
[Q,R] = qr(X);
indx = find(diag(R)<0);
if ~isempty(indx)
    Q(:,indx) = -Q(:,indx);
    R(indx,:) = -R(indx,:);
end