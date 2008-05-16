function dynareroot = dynare_config(path_to_dynare)
% This function tests the existence of valid mex files (for qz
% decomposition, solution to sylvester equation and kronecker
% products...) and, if needed, add paths to the matlab versions
% of these routines.
% 
%
% INPUTS
%   none
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none
%  
% part of DYNARE, copyright Dynare Team (2001-2008)
% Gnu Public License.
MATLAB  = ver('matlab');
%% FIXME:
%% It's not satisfactory to convert string versions into numbers, and to
%% compare these numbers:
%% - conversion will fail if version = 1.2.3
%% - it will give 7.10 < 7.9
VERSION = str2num(MATLAB.Version);
if nargin
    addpath(path_to_dynare);
end
dynareroot = strrep(which('dynare.m'),'dynare.m','');
if (VERSION <= 7.4)
    addpath([dynareroot '../mex/2007a/'])
else
    addpath([dynareroot '../mex/2007b/'])
end
%% Set mex routine names
mex_status = cell(1,3);
mex_status(1,1) = {'mjdgges'};
mex_status(1,2) = {'qz'};
mex_status(1,3) = {'Generalized QZ'};
mex_status(2,1) = {'gensylv'};
mex_status(2,2) = {'gensylv'};
mex_status(2,3) = {'Sylvester equation solution'};
mex_status(3,1) = {'A_times_B_kronecker_C'};
mex_status(3,2) = {'kronecker'};
mex_status(3,3) = {'Kronecker products'};
mex_status(4,1) = {'sparse_hessian_times_B_kronecker_C'};
mex_status(4,2) = {'kronecker'};
mex_status(4,3) = {'Sparse kronecker products'};
number_of_mex_files = size(mex_status,1);
%% Remove some directories from matlab's path. This is necessary if the user has
%% added dynare_v4/matlab with the subfolders. Matlab has to ignore these
%% subfolders if valid mex files exist.
matlab_path = path;
for i=1:number_of_mex_files
    test = strfind(matlab_path,[dynareroot mex_status{i,2}]);
    action = length(test);
    if action
        rmpath([dynareroot mex_status{i,2}]);
        matlab_path = path;
    end
end
%% Test if valid mex files are available, if a mex file is not available
%% a matlab version of the routine is included in the path.
disp(' ')
disp('Configuring Dynare ...')
for i=1:number_of_mex_files
    test = (exist(mex_status{i,1}) == 3);
    if ~test
        addpath([dynareroot mex_status{i,2}]);
        message = '[m]   ';
    else
        message = '[mex] '; 
    end
    disp([ message mex_status{i,3} '.' ])
    disp([])
end
disp(' ')