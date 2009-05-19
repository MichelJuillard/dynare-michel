function test_qmc_sequence()

% Copyright (C) 2009 Dynare Team
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

transform = uint32(0);
dimension = uint32(2);
flag = uint32(3);
seed = uint32(1477);

count = 0;

for i=1:20

    disp('Step 1...')
    qmc_init = qmc_sequence(dimension,flag,seed,transform);
    disp('Done!')

    disp('Step 2...')
    number_of_simulations = uint32(20);
    [Q, qmc_new] = qmc_sequence(qmc_init,number_of_simulations);
    disp('Done!')

    number_of_simulations = uint32(10);

    disp('Step 3...')
    [Q1, qmc_new1] = qmc_sequence(qmc_init,number_of_simulations);
    disp('Done!')

    disp('Step 4...')
    [Q2, qmc_new2] = qmc_sequence(qmc_new1,number_of_simulations);
    disp('Done!')
    
    Q3 = Q - [Q1;Q2] ;
    
    test = max(max(abs(Q3)));
    
    if test<2*eps
        count = count+1;
    end
    
    if i==5 | i== 10
        clear mex;
        close all;
    end

end

disp('qmc_sequence:: No error!')
clear all