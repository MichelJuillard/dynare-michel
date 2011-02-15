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

configstr = ['./configure --with-boost=' [cd '/src/boost_1_45_0'] ...
    ' --with-matlab=' matlabroot ' MATLAB_VERSION=' version('-release')];

cd src/mex/build/matlab;

system(configstr, '-echo');
system('make', '-echo');

cd ../../../..

!cp src/mex/matlab/* mex/matlab/
clear all