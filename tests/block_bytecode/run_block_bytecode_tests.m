## Copyright (C) 2010 Dynare Team
##
## This file is part of Dynare.
##
## Dynare is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Dynare is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

addpath(argv(){1})

if !strcmp(dynare_version(), argv(){2})
  error("Incorrect version of Dynare is being tested")
endif

## Ask gnuplot to create graphics in text mode
## Note that setenv() was introduced in Octave 3.0.2, for compatibility
## with MATLAB
putenv("GNUTERM", "dumb")

for block = 0:1
  for bytecode = 0:1
    ## Recall that solve_algo=7 and stack_solve_algo=2 are not supported
    ## under Octave
    default_solve_algo = 2;
    default_stack_solve_algo = 0;
    if !block && !bytecode
      solve_algos = 0:4;
      stack_solve_algos = 0;
    elseif block && !bytecode
      solve_algos = [0:4 6 8];
      stack_solve_algos = [0 1 3 4];
    else
      solve_algos = [0:6 8];
      stack_solve_algos = [0 1 3:5];
    endif

    for i = 1:length(solve_algos)
      save ws
      if !block && !bytecode && (i == 1)
          run_ls2003(block, bytecode, solve_algos(i), default_stack_solve_algo)
          y_ref = oo_.endo_simul;
          save('test.mat','y_ref');
      else
          run_ls2003(block, bytecode, solve_algos(i), default_stack_solve_algo)
          load('test.mat','y_ref');
          diff = oo_.endo_simul - y_ref;
          assert(abs(diff) <= options_.dynatol);
          endif
      load ws
    endfor
    for i = 1:length(stack_solve_algos)
      save ws
      run_ls2003(block, bytecode, default_solve_algo, stack_solve_algos(i))
      load ws
    endfor
  endfor
endfor

## Local variables:
## mode: Octave
## End:
