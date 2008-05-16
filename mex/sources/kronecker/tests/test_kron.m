function test_kron(test)
    
    
    if ~nargin
        test = 3;
    end
    
    
    if test == 1
    
        percentage_of_non_zero_elements = 10e-4;
        NumberOfVariables = 100;
        NumberOfEquations = 50 ;
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
        
        disp('')
        disp('Computation of A*kron(B,B) with the mex file (v1):')
        tic 
            D1 = sparse_hessian_times_B_kronecker_C(A,B);
        toc
        
        disp('')
        disp('Computation of A*kron(B,B) with the mex file (v2):')
        tic 
            D2 = sparse_hessian_times_B_kronecker_C(A,B,B);
        toc
        
        
        disp('');
        disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
        
        
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


        disp(' ')
        disp('Direct computation of A*kron(B,B):')
        tic
            try
                D4 = A*kron(B,B);
                notest = 0;
            catch
                notest = 1;
                disp('Out of memory')
            end
        toc
        if ~notest
            disp('');
            disp(['Difference between D1 and D4 = ' num2str(max(max(abs(D1-D4))))]);        
        end
end





if test > 1

    hessian = 0;
    load nash_matrices;
    
    r2 = size(zx,1);
    c2 = size(zx,2);
    r1 = size(hessian,1);
    c1 = size(hessian,2);
    
    disp(['Percentage of non zero elements in the hessian matrix = ' num2str(100*nnz(hessian)/(c1*c2)) '%'])
    
    disp('')
    disp('Computation of A*kron(B,B) with the mex file (v1):')
    tic 
        D1 = sparse_hessian_times_B_kronecker_C(hessian,zx);
    toc
 
    disp('')
    disp('Computation of A*kron(B,B) with the mex file (v2):')
    tic 
        D2 = sparse_hessian_times_B_kronecker_C(hessian,zx,zx);
    toc
    
    disp('');
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    
    disp(' ')
% $$$     disp('Computation of A*kron(B,B) with two nested loops:')
% $$$     tic
% $$$         D3 = zeros(r1,c2*c2);
% $$$         k = 0;
% $$$         for i1 = 1:c2
% $$$             for i2 = 1:c2
% $$$                 k = k+1;
% $$$                 D3(:,k) = hessian*kron(zx(:,i1),zx(:,i2)); 
% $$$             end
% $$$         end
% $$$     toc
% $$$     
% $$$     disp(' ')
% $$$     disp(['Maximum absolute difference = ' num2str(max(max(abs(D1-D3))))])

    disp(' ')
    disp(['Percentage of non zero elements in the result matrix = ' num2str(100*nnz(D1)/(r1*c2^2)) '%']);
    
end

if test>2
    A = randn(100,100);
    B = randn(10,10);
    C = randn(10,10);
    disp('Test with full format matrix -- 1')
    D1 = A*kron(B,C);
    D2 = A_times_B_kronecker_C(A,B,C);
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
    disp('Test with full format matrix -- 1')
    D1 = A*kron(B,B);
    D2 = A_times_B_kronecker_C(A,B);
    disp(['Difference between D1 and D2 = ' num2str(max(max(abs(D1-D2))))]);
end