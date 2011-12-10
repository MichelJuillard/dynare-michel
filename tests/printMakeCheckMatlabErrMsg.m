function printMakeCheckMatlabErrMsg(modfilename, exception)
    fprintf('\n********************************************\n');
    disp('*** DYNARE-TEST-MATLAB ERROR ENCOUNTERED ***');
    disp('********************************************');
    disp(['  WHILE RUNNING MODFILE: ' modfilename]);
    fprintf('\n');
    disp(getReport(exception));
    fprintf('*************************************\n\n\n');
end
