% Copyright (C) 2011-2013 Dynare Team
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

% Implementation notes:
%
% Before every call to Dynare, the contents of the workspace is saved in
% 'wsMat.mat', and reloaded after Dynare has finished (this is necessary since
% Dynare does a 'clear -all').  Also note that we take care of clearing the
% 'exception' variable in all 'catch' block, because otherwise the next 'load
% wsMat' within a 'catch' block will overwrite the last exception.

top_test_dir = getenv('TOP_TEST_DIR');
addpath(top_test_dir);
addpath([top_test_dir filesep '..' filesep 'matlab']);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% Test block_bytecode/ls2003.mod with various combinations of
% block/bytecode/solve_algo/stack_solve_algo
failedBlock = {};
num_block_tests = 0;
cd([top_test_dir filesep 'block_bytecode']);
has_optimization_toolbox = user_has_matlab_license('optimization_toolbox');
for blockFlag = 0:1
    for bytecodeFlag = 0:1
        default_solve_algo = 2;
        default_stack_solve_algo = 0;
        if ~blockFlag && ~bytecodeFlag
            solve_algos = 1:4;
            stack_solve_algos = [0 6];
        elseif blockFlag && ~bytecodeFlag
            solve_algos = [1:4 6:8];
            stack_solve_algos = 0:4;
        else
            solve_algos = 1:8;
            stack_solve_algos = 0:5;
        end
        if has_optimization_toolbox
            solve_algos = [ solve_algos 0 ];
        end

        for i = 1:length(solve_algos)
            num_block_tests = num_block_tests + 1;
            if ~blockFlag && ~bytecodeFlag && (i == 1)
                % This is the reference simulation path against which all
                % other simulations will be tested
                try
                    old_path = path;
                    save wsMat
                    run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
                    load wsMat
                    path(old_path);
                    y_ref = oo_.endo_simul;
                    save('test.mat','y_ref');
                catch exception
                    load wsMat
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                    clear exception
                end
            else
                try
                    old_path = path;
                    save wsMat
                    run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
                    load wsMat
                    path(old_path);
                    % Test against the reference simulation path
                    load('test.mat','y_ref');
                    diff = oo_.endo_simul - y_ref;
                    if(abs(diff) > options_.dynatol.x)
                        failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                        exception = MException('ERROR: simulation path differs from the reference path');
                        printMakeCheckMatlabErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                        clear exception
                    end
                catch exception
                    load wsMat
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                    clear exception
                end
            end
        end
        for i = 1:length(stack_solve_algos)
            num_block_tests = num_block_tests + 1;
            try
                old_path = path;
                save wsMat
                if blockFlag && ~bytecodeFlag && stack_solve_algos(i) == 3
                    error('This test currently enters an infinite loop, skipping')
                end
                run_ls2003(blockFlag, bytecodeFlag, default_solve_algo, stack_solve_algos(i))
                load wsMat
                path(old_path);
                % Test against the reference simulation path
                load('test.mat','y_ref');
                diff = oo_.endo_simul - y_ref;
                if(abs(diff) > options_.dynatol.x)
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                    exception = MException('ERROR: simulation path difers from the reference path');
                    printMakeCheckMatlabErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], exception);
                    clear exception
                end
            catch exception
                load wsMat
                path(old_path);
                failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                printMakeCheckMatlabErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(bytecodeFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], exception);
                clear exception
            end
        end
    end
end
delete('wsMat.mat')
cd(getenv('TOP_TEST_DIR'));
fid = fopen('run_block_byte_tests_matlab.m.trs', 'w+');
if size(failedBlock,2) > 0
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: %d\n', num_block_tests);
  fprintf(fid,':number-failed-tests: %d\n', size(failedBlock,2));
  fprintf(fid,':list-of-failed-tests: %s\n', failedBlock{:});
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: %d\n', num_block_tests);
  fprintf(fid,':number-failed-tests: 0\n');
end
fclose(fid);
exit;
