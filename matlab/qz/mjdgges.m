function [ss,tt,w,sdim,eigval,info] = mjdgges(e,d,qz_criterium)
% QZ decomposition, Sims' codes are used.
%
% INPUTS
%   e            [double] real square (n*n) matrix.
%   d            [double] real square (n*n) matrix.
%   qz_criterium [double] scalar (1+epsilon).
%    
% OUTPUTS
%   ss           [complex] (n*n) matrix.
%   tt           [complex] (n*n) matrix.
%   w            [complex] (n*n) matrix.
%   sdim         [integer] scalar.    
%   eigval       [complex] (n*1) vector. 
%   info         [integer] scalar.
%    
% ALGORITHM
%   Sims's qzdiv routine is used.
%
% SPECIAL REQUIREMENTS
%   none.
%  
%  
% part of DYNARE, copyright Dynare Team (1996-2008)
% Gnu Public License.
    
% Chek number of inputs and outputs.
if nargin>3 | nargin<2
    error('Three or two input arguments required!')
end
if nargout>6
    error('No more than six output arguments!')
end
% Check the first two inputs.
[me,ne] = size(e);
[md,nd] = size(d);
if ( ~isreal(e) | ~isreal(d) | me~=ne | md~=nd | me~=nd)
    % info should be negative in this case, see dgges.f.
    error('MYDGGES requires two square real matrices of the same dimension.')
end
% Set default value of qz_criterium.
if nargin <3
    qz_criterium = 1 + 1e-6; 
end
% Initialization of the output arguments.
ss = zeros(ne,ne);
tt = zeros(ne,ne);
w  = zeros(ne,ne);
sdim   = 0;
eigval = zeros(ne,1);
info   = 0;
% Computational part.
try
    [ss,tt,qq,w] = qz(e,d);
    [tt,ss,qq,w] = qzdiv(qz_criterium,tt,ss,qq,w);
    warning_old_state = warning;
    warning off;
    eigval = diag(ss)./diag(tt);
    warning warning_old_state;
    sdim = sum( abs(eigval) < qz_criterium );
catch
    info = 1;% Not as precise as lapack's info!
end