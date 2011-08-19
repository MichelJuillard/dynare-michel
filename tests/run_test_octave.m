## Copyright (C) 2009-2011 Dynare Team
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

## Test Dynare Version
if !strcmp(dynare_version(), dynver)
  error("Incorrect version of Dynare is being tested")
endif

## List of files to be tested
name = filesToTest();

## Ask gnuplot to create graphics in text mode
## Note that setenv() was introduced in Octave 3.0.2, for compatibility
## with MATLAB
putenv("GNUTERM", "dumb")

## BASE TESTS
failedBase = {};

top_test_dir = pwd;
addpath(top_test_dir);
for i=1:size(name,2)
  try
    [directory, testfile, ext] = fileparts([top_test_dir '/' name{i}]);
    cd(directory);
    printf("***  TESTING: %s ***\n", name{i});
    dynare([testfile ext], 'noclearall')
  catch
    failedBase{size(failedBase,2)+1} = name{i};
    printMakeCheckOctaveErrMsg(name{i}, lasterror);
  end_try_catch

  cd(top_test_dir);
  save('makeCheckOctaveBase.mat', 'name', 'failedBase', 'i', 'top_test_dir');
  clear -all;
  load('makeCheckOctaveBase.mat');
end

## BLOCK TEST
clear i name;
failedBlock = {};
num_block_tests = 0;
cd([top_test_dir '/block_bytecode']);
save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
for blockFlag = 0:1
  for bytecodeFlag = 0:1
    ## Recall that solve_algo=7 and stack_solve_algo=2 are not supported
    ## under Octave
    default_solve_algo = 2;
    default_stack_solve_algo = 0;
    if !blockFlag && !bytecodeFlag
      solve_algos = 0:4;
      stack_solve_algos = 0;
    elseif blockFlag && !bytecodeFlag
      solve_algos = [0:4 6 8];
      stack_solve_algos = [0 1 3 4];
    else
      solve_algos = [0:6 8];
      stack_solve_algos = [0 1 3:5];
    endif

    for i = 1:length(solve_algos)
      num_block_tests = num_block_tests + 1;
      save wsOct
      if !blockFlag && !bytecodeFlag && (i == 1)
        try
          run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
          y_ref = oo_.endo_simul;
          save('test.mat','y_ref');
        catch
          load wsOct
          load('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
          failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
          printMakeCheckOctaveErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
          save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
        end_try_catch
      else
        try
          run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
          load('test.mat','y_ref');
          diff = oo_.endo_simul - y_ref;
          if(abs(diff) > options_.dynatol)
            load wsOct
            load('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
            failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
            differr.message = ["makecheck found error: if (abs(diff) <= options_.dynatol) FAILS." ];
            differr.stack(1).file = "run_test_octave.m";
            differr.stack(1).name = "run_test_octave.m";
            differr.stack(1).line = 96;
            differr.stack(1).column = 1;
            printMakeCheckOctaveErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], differr);
            save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
          endif
        catch
          load wsOct
          load('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
          failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
          printMakeCheckOctaveErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
          save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
        end_try_catch
      endif
      load wsOct
    endfor
    for i = 1:length(stack_solve_algos)
      num_block_tests = num_block_tests + 1;
      save wsOct
      try
        run_ls2003(blockFlag, bytecodeFlag, default_solve_algo, stack_solve_algos(i))
      catch
        load wsOct
        load('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
        failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
        printMakeCheckOctaveErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
        save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir');
      end_try_catch
      load wsOct
    endfor
  endfor
endfor

load('makeCheckOctaveBlockByte.mat');
save('makeCheckOctaveBlockByte.mat', 'failedBlock', 'top_test_dir', 'num_block_tests');
delete('wsOct');
clear -all;

load('makeCheckOctaveBlockByte.mat');
delete('makeCheckOctaveBlockByte.mat');

cd(top_test_dir);
load('makeCheckOctaveBase.mat');
delete('makeCheckOctaveBase.mat');

total_tests = size(name,2)+num_block_tests;

% print output to screen and to file
fid = fopen("run_test_octave_output.txt", "w");

printf("\n\n\n");
fprintf(fid,'\n\n\n');
printf("***************************************\n");
fprintf(fid,"***************************************\n");
printf("*         DYNARE TEST RESULTS         *\n");
fprintf(fid,"*         DYNARE TEST RESULTS         *\n");
printf("*        for make check-octave        *\n");
fprintf(fid,"*        for make check-octave        *\n");
printf("***************************************\n");
fprintf(fid,"***************************************\n");
printf("  %d tests PASSED out of %d tests run\n", total_tests-size(failedBase,2)-size(failedBlock,2), total_tests);
fprintf(fid," %d tests PASSED out of %d tests run\n", total_tests-size(failedBase,2)-size(failedBlock,2), total_tests);
printf("***************************************\n");
fprintf(fid,"***************************************\n");
if size(failedBase,2) > 0 || size(failedBlock,2) > 0
  printf("List of %d tests FAILED:\n", size(failedBase,2)+size(failedBlock,2));
  fprintf(fid,"List of %d tests FAILED:\n", size(failedBase,2)+size(failedBlock,2));
  for i=1:size(failedBase,2)
    printf("   * %s\n",failedBase{i});
    fprintf(fid,"   * %s\n", failedBase{i});
  end
  for i=1:size(failedBlock,2)
    printf("   * %s\n",failedBlock{i});
    fprintf(fid,"   * %s\n", failedBlock{i});
  end
  printf("***************************************\n\n");
  fprintf(fid,"***************************************\n\n");
  clear -all
  error("make check-octave FAILED");
end
fclose(fid);
clear -all

## Local variables:
## mode: Octave
## End:
