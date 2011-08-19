% Copyright (C) 2011 Dynare Team
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

% List of files to be tested
name=filesToTest();

% BASE TESTS
failedBase = {};

top_test_dir = pwd;
addpath(top_test_dir);
addpath([top_test_dir '/../matlab']);
for i=1:size(name,2)
    try
        [directory, testfile, ext] = fileparts([top_test_dir '/' name{i}]);
        cd(directory);
        disp(['***  TESTING: ' name{i} ' ***']);
        dynare([testfile ext], 'noclearall')
    catch exception
        failedBase{size(failedBase,2)+1} = name{i};
        printMakeCheckMatlabErrMsg(name{i}, exception);
    end

    cd(top_test_dir);
    save('makeCheckMatlabBase.mat', 'name', 'failedBase', 'i', 'top_test_dir');
    clear -all;
    load('makeCheckMatlabBase.mat');
end

% BLOCK TEST
clear i name;
failedBlock = {};
num_block_tests = 0;
cd([top_test_dir '/block_bytecode']);
save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
for blockFlag = 0:1
    for bytecodeFlag = 0:1
        default_solve_algo = 2;
        default_stack_solve_algo = 0;
        if ~blockFlag && ~bytecodeFlag
            solve_algos = 0:4;
            stack_solve_algos = 0;
        elseif blockFlag && ~bytecodeFlag
            solve_algos = [0:4 6:8];
            stack_solve_algos = 0:4;
        else
            solve_algos = 0:8;
            stack_solve_algos = 0:5;
        end

        for i = 1:length(solve_algos)
            num_block_tests = num_block_tests + 1;
            save wsMat
            if ~blockFlag && ~bytecodeFlag && (i == 1)
                try
                    run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
                    y_ref = oo_.endo_simul;
                    save('test.mat','y_ref');
                catch exception
                    load wsMat
                    load('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                    save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                end
            else
                try
                    run_ls2003(blockFlag, bytecodeFlag, solve_algos(i), default_stack_solve_algo)
                    load('test.mat','y_ref');
                    diff = oo_.endo_simul - y_ref;
                    if(abs(diff) > options_.dynatol)
                        load wsMat
                        load('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                        failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                        exception = MException('ERROR: diff fails for block code in run_test_octave.m', ['makecheck found error: if (abs(diff) <= options_.dynatol) FAILS. line 85, col 1.' ]);
                        printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                        save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                    end
                catch exception
                    load wsMat
                    load('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                    save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                end
            end
            load wsMat
        end
        for i = 1:length(stack_solve_algos)
            num_block_tests = num_block_tests + 1;
            save wsMat
            try
                run_ls2003(blockFlag, bytecodeFlag, default_solve_algo, stack_solve_algos(i))
            catch exception
                load wsMat
                load('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
                failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir');
            end
            load wsMat
        end
    end
end

load('makeCheckBlockByteMatlab.mat');
save('makeCheckBlockByteMatlab.mat', 'failedBlock', 'top_test_dir', 'num_block_tests');
delete('wsMat.mat');
clear -all;

load('makeCheckBlockByteMatlab.mat');
delete('makeCheckBlockByteMatlab.mat');

cd(top_test_dir);
load('makeCheckMatlabBase.mat');
delete('makeCheckMatlabBase.mat');

total_tests = size(name,2)+num_block_tests;

% print output to screen and to file
fid = fopen('run_test_matlab_output.txt', 'w');

fprintf('\n\n\n');
fprintf(fid,'\n\n\n');
disp('***************************************');
fprintf(fid,'***************************************\n');
disp('*         DYNARE TEST RESULTS         *');
fprintf(fid,'*         DYNARE TEST RESULTS         *\n');
disp('*        for make check-matlab        *');
fprintf(fid,'*        for make check-matlab        *\n');
disp('***************************************');
fprintf(fid,'***************************************\n');
disp(['  ' num2str(total_tests-size(failedBase,2)-size(failedBlock,2)) ' tests PASSED out of ' num2str(total_tests) ' tests run']);
fprintf(fid,' %d tests PASSED out of %d tests run\n', total_tests-size(failedBase,2)-size(failedBlock,2), total_tests);
disp('***************************************');
fprintf(fid,'***************************************\n');
if size(failedBase,2) > 0 || size(failedBlock,2) > 0
    disp(['List of ' num2str(size(failedBase,2)+size(failedBlock,2)) ' tests FAILED:']);
    fprintf(fid,'List of %d tests FAILED:\n', size(failedBase,2)+size(failedBlock,2));
    for i=1:size(failedBase,2)
        disp(['   * ' failedBase{i}]);
        fprintf(fid,'   * %s\n', failedBase{i});
    end
    for i=1:size(failedBlock,2)
        disp(['   * ' failedBlock{i}]);
        fprintf(fid,'   * %s\n', failedBlock{i});
    end
    fprintf('***************************************\n\n');
    fprintf(fid,'***************************************\n\n');
end
fclose(fid);
exit;
