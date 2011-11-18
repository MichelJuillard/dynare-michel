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

% Implementation notes:
%
% Before every call to Dynare, the contents of the workspace is saved in
% 'wsMat.mat', and reloaded after Dynare has finished (this is necessary since
% Dynare does a 'clear -all').  Also note that we take care of clearing the
% 'exception' variable in all 'catch' block, because otherwise the next 'load
% wsMat' within a 'catch' block will overwrite the last exception.

top_test_dir = pwd;
addpath(top_test_dir);
addpath([top_test_dir '/../matlab']);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% Test MOD files listed in Makefile.am
name=getenv('MODFILES');
num_modfiles = 0;

failedBase = {};

while ~isempty(name)
    [modfile, name] = strtok(name);
    num_modfiles = num_modfiles + 1;
    [directory, testfile, ext] = fileparts([top_test_dir '/' modfile]);
    cd(directory);
    disp('');
    disp(['***  TESTING: ' modfile ' ***']);
    try
        old_path = path;
        save wsMat
        dynare([testfile ext],'console')
        clear -all
        load wsMat
        path(old_path);
    catch exception
        clear -all
        load wsMat
        path(old_path);
        failedBase{size(failedBase,2)+1} = modfile;
        printMakeCheckMatlabErrMsg(modfile, exception);
        clear exception
    end
    close all
    delete('wsMat.mat')

    cd(top_test_dir);
end

% Test block_bytecode/ls2003.mod with various combinations of
% block/bytecode/solve_algo/stack_solve_algo
failedBlock = {};
num_block_tests = 0;
cd([top_test_dir '/block_bytecode']);
for blockFlag = 0:1
    for bytecodeFlag = 0:1
        default_solve_algo = 2;
        default_stack_solve_algo = 0;
        if ~blockFlag && ~bytecodeFlag
            solve_algos = 1:4;
            stack_solve_algos = 0;
        elseif blockFlag && ~bytecodeFlag
            solve_algos = [1:4 6:8];
            stack_solve_algos = 0:4;
        else
            solve_algos = 1:8;
            stack_solve_algos = 0:5;
        end
        if license('test', 'optimization_toolbox')
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
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
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
                    if(abs(diff) > options_.dynatol)
                        failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                        exception = MException('ERROR: simulation path differs from the reference path');
                        printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                        clear exception
                    end
                catch exception
                    load wsMat
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                    clear exception
                end
            end
        end
        for i = 1:length(stack_solve_algos)
            num_block_tests = num_block_tests + 1;
            try
                old_path = path;
                save wsMat
                run_ls2003(blockFlag, bytecodeFlag, default_solve_algo, stack_solve_algos(i))
                load wsMat
                path(old_path);
                % Test against the reference simulation path
                load('test.mat','y_ref');
                diff = oo_.endo_simul - y_ref;
                if(abs(diff) > options_.dynatol)
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(default_solve_algo) ', ' num2str(stack_solve_algos(i)) ')'];
                    exception = MException('ERROR: simulation path difers from the reference path');
                    printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(default_solve_algo) ', ' num2str(stack_solve_algos(i)) ')'], exception);
                    clear exception
                end
            catch exception
                load wsMat
                path(old_path);
                failedBlock{size(failedBlock,2)+1} = ['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'];
                printMakeCheckMatlabErrMsg(['block_bytecode/run_ls2003.m(' num2str(blockFlag) ', ' num2str(bytecodeFlag) ', ' num2str(solve_algos(i)) ', ' num2str(default_stack_solve_algo) ')'], exception);
                clear exception
            end
        end
    end
end
delete('wsMat.mat')

cd(top_test_dir);

total_tests = num_modfiles+num_block_tests;

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
