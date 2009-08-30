function b = fn_tran_a2b(A0,Ui,nvar,n0)
% b = fn_tran_a2b(A0,Ui,nvar,n0)
%    Transform A0 to free parameters b's.  Note: columns correspond to equations
%    See Waggoner and Zha's ``A Gibbs sampler for structural VARs''
%
% A0: nvar-by-nvar, contempareous matrix (columns correspond to equations)
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
% nvar:  number of endogeous variables
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
%----------------
% b: sum(n0)-by-1 vector of A0 free parameters
%
% Tao Zha, February 2000.  Revised, August 2000

n0=n0(:);
n0cum = [0; cumsum(n0)];
b=zeros(n0cum(end),1);
for kj = 1:nvar
   b(n0cum(kj)+1:n0cum(kj+1))=Ui{kj}'*A0(:,kj);
end
