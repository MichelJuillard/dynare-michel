function info = test_kron(test,number_of_threads)
% Copyright (C) 2007-2010 Dynare Team
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

if ~nargin
    test = 2;
end

info = 1;

if test == 1
    test1_1 = 1;
    test1_2 = 1;
    disp('')
    disp('I''m building the test problem...')
    tic
    percentage_of_non_zero_elements = 10e-4;
    NumberOfVariables = 549;%100;
    NumberOfEquations = 256;%50
    NumberOfColsInB   = 50 ; 
    A = zeros(NumberOfEquations,NumberOfVariables^2);
    for i = 1:NumberOfEquations
        for j = 1:NumberOfVariables
            for k = j:NumberOfVariables
                if rand<percentage_of_non_zero_elements
                    A(i,(j-1)*NumberOfVariables+k) = randn;
                end
            end
            for h = j+1:NumberOfVariables
                A(i,(h-1)*NumberOfVariables+j) = A(i,(j-1)*NumberOfVariables+h); 
            end
        end
    end
    A = sparse(A);
    B = randn(NumberOfVariables,NumberOfColsInB);
    disp('Done!')
    toc
    disp('')
    disp('Computation of A*kron(B,B) with the mex file (v1):')
    tic 
    [err, D1] = sparse_hessian_times_B_kronecker_C(A,B,number_of_threads);
    mexErrCheck('sparse_hessian_times_B_kronecker_C', err);
    toc
    disp('')
    disp('Computation of A*kron(B,B) with the mex file (v2):')
    tic
    [err, D2] = sparse_hessian_times_B_kronecker_C(A,B,B,number_of_threads);
    mexErrCheck('sparse_hessian_times_B_kronecker_C', err);
    toc
    disp('');
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    if max(max(abs(D1-D2)))>1e-10
        test1_1=0; 
    end
    disp(' ')
    disp('Computation of A*kron(B,B) with two nested loops:')
    tic
    D3 = zeros(NumberOfEquations,NumberOfColsInB*NumberOfColsInB);
    k = 1;
    for i1 = 1:NumberOfColsInB
        for i2 = 1:NumberOfColsInB
            D3(:,k) = A*kron(B(:,i1),B(:,i2)); 
            k = k+1;
        end
    end
    toc
    disp('');
    disp(['Difference between D1 and D3 = ' num2str(max(max(abs(D1-D3))))]);
    if max(max(abs(D1-D3)))>1e-10
        test1_2=0; 
    end
% $$$ FOR THE DIMENSIONS CONSIDERED HERE THIS PART WILL RESULT IN A OUT OF MEMORY ERROR.   
% $$$     disp(' ')
% $$$     disp('Direct computation of A*kron(B,B):')
% $$$     tic
% $$$     try
% $$$         D4 = A*kron(B,B);
% $$$         notest = 0;
% $$$     catch
% $$$         notest = 1;
% $$$         disp('Out of memory')
% $$$     end
% $$$     toc
% $$$     if ~notest
% $$$         disp('');
% $$$         disp(['Difference between D1 and D4 = ' num2str(max(max(abs(D1-D4))))]);        
% $$$     end
    if ~(test1_1 && test1_2)
        info = 0;
    end
end

if test==2
    test2_1 = 1;
    test2_2 = 1;
    disp('I''m loading a real life problem...')
    hessian = 0;
    load nash_matrices;
    disp('Done!')
    r2 = size(zx,1);
    c2 = size(zx,2);
    r1 = size(hessian,1);
    c1 = size(hessian,2);
    disp(' ')
    disp(['Percentage of non zero elements in the hessian matrix = ' num2str(100*nnz(hessian)/(c1*c2)) '%'])
    disp(' ')
    disp('Computation of A*kron(B,B) with the mex file (v1):')
    tic
    [err, D1] = sparse_hessian_times_B_kronecker_C(hessian,zx,number_of_threads);
    mexErrCheck('sparse_hessian_times_B_kronecker_C', err);
    toc
    disp(' ')
    disp('Computation of A*kron(B,B) with the mex file (v2):')
    tic
    [err, D2] = sparse_hessian_times_B_kronecker_C(hessian,zx,zx,number_of_threads);
    mexErrCheck('sparse_hessian_times_B_kronecker_C', err);    
    toc
    disp(' ');
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    if max(max(abs(D1-D2)))>1e-10
        test2_1=0; 
    end
    disp(' ')
    disp('Computation of A*kron(B,B) with two nested loops:')
    tic
    D3 = zeros(r1,c2*c2);
    k = 0;
    for i1 = 1:c2
        for i2 = 1:c2
            k = k+1;
            D3(:,k) = hessian*kron(zx(:,i1),zx(:,i2)); 
        end
    end
    toc
    disp(' ')
    disp(['Maximum absolute difference = ' num2str(max(max(abs(D1-D3))))])
    if max(max(abs(D1-D3)))>1e-10
        test2_2=0; 
    end
    disp(' ')
    disp(['Percentage of non zero elements in the result matrix = ' num2str(100*nnz(D1)/(r1*c2^2)) '%']);
    if ~(test2_1 && test2_2)
        info = 0;
    end
end

if test==3
    test3_1 = 1;
    test3_2 = 1;
    A = randn(100,10000);
    B = randn(100,10);
    C = randn(100,10);
    disp('Test with full format matrix -- 1(a)')
    D1 = A*kron(B,C);
    tic
    [err, D2] = A_times_B_kronecker_C(A,B,C,number_of_threads);
    mexErrCheck('A_times_B_kronecker_C', err);
    toc
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    if max(max(abs(D1-D2)))>1e-10
        test3_1=0; 
    end
    disp('Test with full format matrix -- 1(b)')
    D1 = A*kron(B,B);
    tic
    [err, D2] = A_times_B_kronecker_C(A,B,number_of_threads);
    mexErrCheck('A_times_B_kronecker_C', err);
    toc
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    if max(max(abs(D1-D2)))>1e-10
        test3_2=0; 
    end
    if ~(test3_1 && test3_2)
        info = 0;
    end
end