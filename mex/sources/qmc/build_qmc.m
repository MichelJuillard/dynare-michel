addpath '../../../matlab';

MATLAB_PATH = matlabroot;

COMPILE_OPTIONS = '';

% mwSize, mwIndex and mwSignedIndex appeared in Matlab 7.3
if matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMWTYPES_NOT_DEFINED' ];
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer)) ...
      && ~matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
end

if matlab_ver_less_than('7.5')
    OUTPUT_DIR = '../../2007a';
else
    OUTPUT_DIR = '../../2007b';
end

disp(' ')
if exist(OUTPUT_DIR,'dir')
    disp('Delete old qmc mex file.')
    delete([OUTPUT_DIR '/qmc_sequence.' mexext]);
else
    whereami = pwd;
    disp(['Create directory ' whereami(1:end-7) OUTPUT_DIR(4:end) '.'])
    mkdir(OUTPUT_DIR);
end
disp(' ')

% Comment next line to suppress compilation debugging info
% COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v '  ];

% Set Optimization and Debug flags
CXXOPTIMFLAGS = ' CXXOPTIMFLAGS=-O3 ';
COPTIMFLAGS = ' COPTIMFLAGS=-O3 ';
CXXDEBUGFLAGS = ' CXXDEBUGFLAGS= ';
CDEBUGFLAGS = ' CDEBUGFLAGS= ';
LDOPTIMFLAGS = ' LDOPTIMFLAGS=-O3 ';
LDDEBUGFLAGS = ' LDDEBUGFLAGS= ';

COMPILE_OPTIONS = [ COMPILE_OPTIONS CDEBUGFLAGS COPTIMFLAGS CXXDEBUGFLAGS CXXOPTIMFLAGS LDDEBUGFLAGS LDOPTIMFLAGS];

FC='FC=/usr/bin/gfortran ' ;
FFLAGS= 'FFLAGS= ' ;
FLIBS = 'FLIBS=-L/usr/lib -lgfortran ' ;
FOPTIMFLAGS = 'FOPTIMFLAGS=-O ' ;

COMPILE_OPTIONS = [ COMPILE_OPTIONS FC FFLAGS FOPTIMFLAGS FLIBS ];

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ] ;

disp('Compiling qmc')
system('gfortran -c -O3 low_discrepancy.f');
eval([ COMPILE_COMMAND ' qmc_sequence.cc low_discrepancy.o' ]) ;