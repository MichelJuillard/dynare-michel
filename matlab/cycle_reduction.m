function [X, info] = cycle_reduction(A0, A1, A2, cvg_tol, ch)

%@info:
%! @deftypefn {Function File} {[@var{X}, @var{info}] =} cycle_reduction (@var{A0},@var{A1},@var{A2},@var{cvg_tol},@var{ch})
%! @anchor{cycle_reduction}
%! @sp 1
%! Solves the quadratic matrix equation A2*X^2 + A1*X + A0 = 0.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A0
%! Square matrix of doubles, n*n.
%! @item A1
%! Square matrix of doubles, n*n.
%! @item A2
%! Square matrix of doubles, n*n.
%! @item cvg_tol
%! Scalar double, tolerance parameter.
%! @item ch
%! Any matlab object, if not empty the solution is checked.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item X
%! Square matrix of doubles, n*n, solution of the matrix equation.
%! @item info
%! Scalar integer, if nonzero the algorithm failed in finding the solution of the matrix equation.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @strong{References:}
%! @sp 1
%! D.A. Bini, G. Latouche, B. Meini (2002), "Solving matrix polynomial equations arising in queueing problems", Linear Algebra and its Applications 340, pp. 222-244
%! @sp 1
%! D.A. Bini, B. Meini (1996), "On the solution of a nonlinear matrix equation arising in queueing problems", SIAM J. Matrix Anal. Appl. 17, pp. 906-926.
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

max_it = 300;
it = 0;
info = 0;
crit = 1+cvg_tol;
A_0 = A1;
A0_0 = A0;
A1_0 = A1;
A2_0 = A2;
n = length(A0);
id = 1:n;

while crit > cvg_tol && it < max_it;
    tmp = [A2_0; A0_0]/A1_0;
    TMP = tmp*A0_0;
    A2_1 = - tmp(id,:)* A2_0;
    A_1 = A_0 - TMP(id,:);
    A1_1 = A1_0 - tmp(n+id,:) * A2_0 - TMP(id,:);
    crit = sum(sum(abs(A0_0)));
    A_0 = A_1;
    A0_0 = -TMP(n+id,:);
    A1_0 = A1_1;
    A2_0 = A2_1;
    it = it + 1;
end

if it==max_it
    disp(['convergence not achieved after ' int2str(it) ' iterations']);
    info = 1;
end

X = -A_0\A0;

if (nargin == 5 && ~isempty(ch) )
    %check the solution
    res = A0 + A1 * X + A2 * X * X;
    if (sum(sum(abs(res))) > cvg_tol)
        disp(['the norm residual of the residu ' num2str(res) ' compare to the tolerance criterion ' num2str(cvg_tol)]);
    end
end