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

function run_ls2003(block, bytecode, solve_algo, stack_solve_algo)
  printf("\nTEST: ls2003 (block=%d, bytecode=%d, solve_algo=%d, stack_solve_algo=%d)...\n", block, bytecode, solve_algo, stack_solve_algo);
  fid = fopen("ls2003_tmp.mod", "w");
  assert(fid > 0);
  fprintf(fid, "@#define block = %d\n@#define bytecode = %d\n@#define solve_algo = %d\n@#define stack_solve_algo = %d\n@#include \"ls2003.mod\"\n", block, bytecode, solve_algo, stack_solve_algo)
  fclose(fid);
  dynare("ls2003_tmp.mod","noclearall")
endfunction

## Local variables:
## mode: Octave
## End:
