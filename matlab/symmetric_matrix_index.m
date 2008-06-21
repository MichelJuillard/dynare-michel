function k = symmetric_matrix_index(i,j,n)
% part of DYNARE, copyright Dynare Team (2007-2008).
% Gnu Public License.    
    k = (i-1)*n+j-i*(i-1)/2;