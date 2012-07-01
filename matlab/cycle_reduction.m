function [X, info] = cycle_reduction(A0, A1, A2, cvg_tol, ch)
% function [X, info] = cycle_reduction(A0,A1,A2,A3, cvg_tolch)
% 
% Solves Polynomial Equation: 
% A0 + A1 * X + A2 * X² = 0
% Using Cyclic Reduction algorithm
% -  D.A. Bini, G. Latouche, B. Meini (2002), "Solving matrix polynomial equations arising in
%    queueing problems", Linear Algebra and its Applications 340 (2002) 225–244
% -  D.A. Bini, B. Meini, On the solution of a nonlinear matrix equation arising in queueing problems,
%    SIAM J. Matrix Anal. Appl. 17 (1996) 906–926.
% =================================================================

% Copyright (C) 2006-2012 Dynare Team
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
  while crit > cvg_tol && it < max_it;
       i_A1_0 = inv(A1_0);
       A2_0_i_A1_0 = A2_0 * i_A1_0;
       A0_0_i_A1_0 = A0_0 * i_A1_0;
       A1_INC = A2_0_i_A1_0 * A0_0;
       A_1 = A_0 - A1_INC;
       A0_1 = - A0_0_i_A1_0 * A0_0;
       A1_1 = A1_0 - A0_0_i_A1_0 * A2_0 - A1_INC;
       A2_1 = - A2_0_i_A1_0 * A2_0;


      crit = sum(sum(abs(A0_0)));

      A_0 = A_1;
      A0_0 = A0_1;
      A1_0 = A1_1;
      A2_0 = A2_1;
      it = it + 1;
      %disp(['it=' int2str(it) ' crit = ' num2str(crit)]);
  end;
  if it==max_it
      disp(['convergence not achieved after ' int2str(it) ' iterations']);
      info = 1;
  end
  X = - inv(A_0) * A0;
  if (nargin == 5 && ~isempty( ch ) == 1 )
      %check the solution
      res = A0 + A1 * X + A2 * X * X;
      if (sum(sum(abs(res))) > cvg_tol)
          disp(['the norm residual of the residu ' num2str(res) ' compare to the tolerance criterion ' num2str(cvg_tol)]);
      end;
  end;