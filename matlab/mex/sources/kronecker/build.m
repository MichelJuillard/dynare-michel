function build()
% This function builds two mex files which compute A*kron(B,C).
%
% Dynare Team, 2003-2007.

    SOURCE1 = 'sparse_hessian_times_B_kronecker_C.cc ' ;
    SOURCE2 = 'A_times_B_kronecker_C.cc ';
    
    MATLAB_PATH = matlabroot;
    
    if strcmpi('GLNX86',computer)%linux (32 bits).
        MATLAB_PATH = [ matlabroot '/bin' ];
        LIB_PATH = [ MATLAB_PATH '/glnx86' ];
        LIB_NAME = '/libmwblas.so';
        COPY_COMMAND = 'cp  *.mexglx ../../dynare_v4/matlab';
        COMPILE_COMMAND = '/mex ';
        COMPILE_OPTIONS = '';
        CLEAN_COMMAND = 'rm *.mexglx';
    elseif strcmpi('GLNXA64',computer)%linux (64 bits).
        MATLAB_PATH = [matlabroot '/bin'];
        LIB_PATH = [MATLAB_PATH '/glnxa64'];
        LIB_NAME = '/libmwblas.so';
        COPY_COMMAND = 'cp  *.mexa64 ../../dynare_v4/matlab/';
        COMPILE_COMMAND = '/mex ';
        COMPILE_OPTIONS = ' -largeArrayDims ';
        CLEAN_COMMAND = 'rm *.mexa64';
    elseif strcmpi('PCWIN',computer)%windows (32 bits).
        MATLAB_PATH = ['"' matlabroot '\bin'];
        LIB_PATH = ['"' matlabroot '\extern\lib\win32\microsoft"'];
        LIB_NAME = '\libmwblas.lib';
        COPY_COMMAND = 'cp  *.mexw32 ../../dynare_v4/matlab/';
        COMPILE_COMMAND = '\mex" ';
        COMPILE_OPTIONS = '-v ';
        CLEAN_COMMAND = 'del *.mexw32';
    end
      
    system([MATLAB_PATH COMPILE_COMMAND COMPILE_OPTIONS SOURCE1]);
    cmd = [MATLAB_PATH COMPILE_COMMAND COMPILE_OPTIONS SOURCE2 LIB_PATH LIB_NAME]
    system(cmd);
    system(COPY_COMMAND);
    system(CLEAN_COMMAND);