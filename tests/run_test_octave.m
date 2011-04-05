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
name = { ...
        'ramst.mod' ...
        'ramst_a.mod' ...
        'example1.mod' ...
        'example2.mod' ...
        'example1_use_dll.mod' ...
        'example1_with_tags.mod' ...
        't_sgu_ex1.mod' ...
        'osr_example.mod' ...
        'optimal_policy/ramsey.mod' ...
        'optimal_policy/mult_elimination_test.mod' ...
        'discretionary_policy/dennis_1.mod' ...
        'ramst_initval_file.mod' ...
        'ramst_normcdf_and_friends.mod' ...
        'example1_varexo_det.mod' ...
        'predetermined_variables.mod' ...
        'fs2000_nonstationary.mod' ...
        'fs2000_ssfile.mod' ...
        'comments.mod' ...
        'histval_sto.mod' ...
        'histval_det.mod' ...
        'expectation.mod' ...
        'steady_state_operator/standard.mod' ...
        'steady_state_operator/use_dll.mod' ...
        'steady_state_operator/block.mod' ...
        'steady_state_operator/bytecode_test.mod' ...
        'block_bytecode/ireland.mod' ...
        'block_bytecode/ramst_normcdf_and_friends.mod' ...
        'k_order_perturbation/fs2000k2a.mod' ...
        'k_order_perturbation/fs2000k2_use_dll.mod' ...
        'k_order_perturbation/fs2000k_1_use_dll.mod' ...
        'k_order_perturbation/fs2000k3_use_dll.mod' ...
        'k_order_perturbation/fs2000k2_m.mod' ...
        'k_order_perturbation/fs2000k_1_m.mod' ...
        'k_order_perturbation/fs2000k3_m.mod' ...
        'partial_information/PItest3aHc0PCLsimModPiYrVarobsAll.mod' ...
        'partial_information/PItest3aHc0PCLsimModPiYrVarobsCNR.mod' ...
        'arima/mod1.mod' ...
        'arima/mod1a.mod' ...
        'arima/mod1b.mod' ...
        'arima/mod1c.mod' ...
        'arima/mod2.mod' ...
        'arima/mod2a.mod' ...
        'arima/mod2b.mod' ...
        'arima/mod2c.mod' ...
        'fs2000/fs2000.mod' ...
        'fs2000/fs2000a.mod' ...
        'fs2000/fs2000c.mod' ...
        'homotopy/homotopy1_test.mod' ...
        'homotopy/homotopy2_test.mod' ...
        'homotopy/homotopy3_test.mod' ...
        'bvar_a_la_sims/bvar_standalone.mod' ...
        'bvar_a_la_sims/bvar_and_dsge.mod' ...
        'AIM/fs2000x10L9_L.mod' ...
        'AIM/fs2000x10L9_L_AIM.mod' ...
        'AIM/fs2000x10_L9_L.mod' ...
        'AIM/fs2000x10_L9_L_AIM.mod' ...
        'AIM/fs2000_b1L1L.mod' ...
        'AIM/fs2000_b1L1L_AIM.mod' ...
        'AIM/ls2003_2L0L.mod' ...
        'AIM/ls2003_2L0L_AIM.mod' ...
        'AIM/ls2003_2L2L.mod' ...
        'AIM/ls2003_2L2L_AIM.mod' ...
        'conditional_variance_decomposition/example1.mod' ...
        'dsge-var/simul_hybrid.mod' ...
        'dsge-var/dsgevar_forward_calibrated_lambda.mod' ...
        'dsge-var/dsgevar_forward_estimated_lambda.mod' ...
        'external_function/example1_1st_and_2nd_deriv_functions_provided.mod' ...
        'external_function/example1_1st_and_2nd_deriv_functions_provided_dll.mod' ...
        'external_function/example1_1st_deriv_function_provided.mod' ...
        'external_function/example1_1st_deriv_function_provided_dll.mod' ...
        'external_function/example1_no_deriv_functions_provided.mod' ...
        'external_function/example1_no_deriv_functions_provided_dll.mod' ...
        'seeds.mod' ...
        'simul/example1.mod' ...
        };

## Ask gnuplot to create graphics in text mode
## Note that setenv() was introduced in Octave 3.0.2, for compatibility
## with MATLAB
putenv("GNUTERM", "dumb")

## BASE TESTS
failed = {};

top_test_dir = pwd;
addpath(top_test_dir);
addpath([top_test_dir '../matlab']);
for i=1:size(name,2)
  try
    [directory, testfile, ext] = fileparts([top_test_dir '/' name{i}]);
    cd(directory);
    printf("***  TESTING: %s ***\n", name{i});
    dynare([testfile ext], 'noclearall')
  catch
    failed{size(failed,2)+1} = name{i};
    printMakeCheckErrMsg(name{i}, lasterror);
  end_try_catch

  cd(top_test_dir);
  save('makeCheckBase.mat', 'name', 'failed', 'i', 'top_test_dir');
  clear -all;
  load('makeCheckBase.mat');
end

## BLOCK TEST
clear i;
clear name;
failed = {};
num_block_tests = 0;
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
      cd([top_test_dir '/block_bytecode']);
      num_block_tests = num_block_tests + 1;
      if !block && !bytecode && (i == 1)
        try
          run_ls2003(block, bytecode, solve_algos(i), default_stack_solve_algo)
          y_ref = oo_.endo_simul;
          save('test.mat','y_ref');
        catch
          failed{size(failed,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
          printMakeCheckErrMsg(['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
        end_try_catch
      else
        try
          run_ls2003(block, bytecode, solve_algos(i), default_stack_solve_algo)
          load('test.mat','y_ref');
          diff = oo_.endo_simul - y_ref;
          if (abs(diff) <= options_.dynatol)
            failed{size(failed,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
            printMakeCheckErrMsg(['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
          end
        catch
          failed{size(failed,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
          printMakeCheckErrMsg(['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
        end_try_catch
      endif

      cd(top_test_dir);
      save('makeCheckBlock.mat', 'failed', 'num_block_tests', 'top_test_dir', 'i', ...
	   'solve_algos', 'block', 'bytecode', 'default_stack_solve_algo', 'default_solve_algo', 'stack_solve_algos');
      clear -all;
      load('makeCheckBlock.mat');
    endfor

    for i = 1:length(stack_solve_algos)
      cd([top_test_dir '/block_bytecode']);
      num_block_tests = num_block_tests + 1;
      try
        run_ls2003(block, bytecode, default_solve_algo, stack_solve_algos(i))
      catch
        failed{size(failed,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
        printMakeCheckErrMsg(['block_bytecode/run_ls2003.m(' num2str(block) ', ' num2str(bytecode) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], lasterror);
      end_try_catch

      cd(top_test_dir);
      save('makeCheckBlock.mat', 'failed', 'num_block_tests', 'top_test_dir', 'i', ...
	   'solve_algos', 'block', 'bytecode', 'default_stack_solve_algo', 'default_solve_algo', 'stack_solve_algos');
      clear -all;
      load('makeCheckBlock.mat');
    endfor
  endfor
endfor

cd(top_test_dir);
clear -all

load('makeCheckBase.mat');
failedBase = failed;
delete('makeCheckBase.mat');

load('makeCheckBlock.mat');
failedBlock = failed;
delete('makeCheckBlock.mat');

total_tests = size(name,2)+num_block_tests;

printf("\n\n\n");
printf("***************************************\n");
printf("*         DYNARE TEST RESULTS         *\n");
printf("***************************************\n");
printf("  %d tests PASSED out of %d tests run\n", total_tests-size(failedBase,2)-size(failedBlock,2), total_tests);
printf("***************************************\n");
if size(failedBase,2) > 0 || size(failedBlock,2) > 0
  printf("List of %d tests FAILED:\n", size(failedBase,2)+size(failedBlock,2));
  for i=1:size(failedBase,2)
    printf("   * %s\n",failedBase{i});
  end
  for i=1:size(failedBlock,2)
    printf("   * %s\n",failedBlock{i});
  end
  printf("***************************************\n\n");
  clear -all
  error("Make Check FAILED");
end
clear -all

## Local variables:
## mode: Octave
## End:
