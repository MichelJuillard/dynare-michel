%
% GENSYLV solves the following matrix equation:
%           A*X + [0 B]*X*kron(C,..,C) = D,
%  where
%       A ...... regular (n,n) matrix,
%       [0 B] .. (n,n) matrix with a few first
%                columns equal to zeros
%       B ...... rectangular (n, nz) matrix with nz<=n 
%       C ...... regular (m,m) matrix, whose spectrum is
%                within (-1, 1)
%       kron(C,..,C)
%         ...... Kronecker power of C of order 'order'
%       D .....  rectangular (n, m^order) matrix.
%
% [err, X] = gensylv(order, A, B, C, D)
%       returns err a scalar where 1 indicates failure, 0 indicates success
%       and X as the solution, doesn't perform any checks
%
% [err, X, par] = gensylv(order, A, B, C, D)
%       solves the system, and performs checks. 'par' is struct
%       containing information about solution and error norms
%       returned by the check. This is a list of the struct
%       members, some of them may be missing in actual returned
%       value:
%       method     = method used for solution recursive/iterative
%       convergence_tol = convergence tolerance for iter. method
%       max_num_iter    = max number of steps for iter. method
%       bs_norm    = Bavely Stewart log10 norm for diagonalization
%       converged  = convergence status for iterative method
%       iter_last_norm  = residual norm of the last step of iterations
%       num_iter   = number of iterations performed
%       f_err1     = norm 1 of abs. error C-V*F*inv(V)
%       f_errI     = norm Inf of abs. error C-V*F*inv(V)
%       viv_err1   = norm 1 of abs. error I-V*inv(V)
%       viv_errI   = norm Inf of abs. error I-V*inv(V)
%       ivv_err1   = norm 1 of abs. error I-inv(V)*V
%       ivv_errI   = norm Inf of abs. error I-inv(V)*V
%       f_blocks   = number of diagonal blocks of F
%       f_largest  = size of largest diagonal block in F
%       f_zeros    = number of off diagonal zeros in F
%       f_offdiag  = number of all offdiagonal elements in F
%       rcondA1    = reciprocal cond 1 number of A
%       rcondAI    = reciprocal cond Inf number of A
%       eig_min    = minimum eigenvalue of vectorized system
%       mat_err1   = rel. matrix 1 norm of A*X-[0 B]*X*kron(C,..,C)-D
%       mat_errI   = rel. matrix Inf norm of       --"--
%       mat_errF   = rel. matrix Frobenius norm of --"--
%       vec_err1   = rel. vector 1 norm of         --"--
%       vec_errI   = rel. vector Inf norm of       --"--
%       cpu_time   = CPU time needed for solution in CPU seconds
%

% $Header: /var/lib/cvs/dynare_cpp/sylv/matlab/gensylv.m,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $
% Tag $Name:  $

function [err, X, varargout] = gensylv(order, A, B, C, D)

% in Windows, ensure that aa_gensylv.dll is loaded, this prevents
% clearing the function and a successive Matlab crash
if strcmp('PCWIN', computer)
  if ~ libisloaded('aa_gensylv') 
    loadlibrary('aa_gensylv', 'dummy');
  end
end

% launch aa_gensylv
if nargout == 2
  X = aa_gensylv(order, A, B, C, D);
elseif nargout == 3
  [X, par] = aa_gensylv(order, A, B, C, D);
  varargout(1) = {par};
else
  error('Must have 2 or 3 output arguments.');
end
err = 0;
  