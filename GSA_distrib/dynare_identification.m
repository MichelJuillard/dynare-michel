In the  as2001 model, after H and J are computed,

rank(H)==size(H,2)
rank(JJ)==size(JJ,2)

vnorm(H)    % vnorm computes the norms of each column of the matrix. I have
attached it below
vnorm(JJ)

% here the problems is due to parameters 9 and 10 which affect only the
steady state

ind1 = find(vnorm(H)~=0);
H1 = H(:,ind1);
JJ1 = JJ(:,ind1);

rank(H1)==size(H1,2)
rank(JJ1)==size(JJ1,2)

% to find near linear dependence problems  I use

for ii = 1:size(H1,2);
    McoH(ii,:) = [ii,cosn([H1(:,ii),H1(:,find([1:1:size(H1,2)]~=ii))])];
    McoJ(ii,:) = [ii,cosn([JJ1(:,ii),JJ1(:,find([1:1:size(JJ1,2)]~=ii))])];
end

format long  % some are nearly 1
McoJ


% here there is no exact linear dependence, but there are several
near-dependencies, mostly due to strong pairwise colliniearities, which can
be checked using

for ii = 1:size(H1,2);
    for jj = 1:size(H1,2);
    PcoH(ii,jj) = [cosn([H1(:,ii),H1(:,jj)])];
    PcoJ(ii,jj) = [cosn([JJ1(:,ii),JJ1(:,jj)])];
    end
end


the cosn approach is more useful for near linear dependence, i.e. weak
identification.  For checking identification only, using the eigenvalues of
(H'H) is easier, because it tells you immediately which columns are
linearly dependent. For instance if 1,3 and 5 columns are linearly
dependent, this can be discovered by  [e1,e2]= eig(H'*H);  first column of
e1, corresponding to the 0 eigenvalue (the (1,1) % element of s2), has
non-zero 1,3 and 5 element.   I think at this point we cannot do much
regarding weak identifiication in dynare, except perhaps some indications
such as small norms of the columns of H or J and strong but not perfect
collinearity.


My suggestion is to have the following steps for identification check in
dynare:

1. check rank of H at theta
      if rank(H)<length(theta) - find out which parameters are involved,
using something like the vnorm and the eigenvalue decomposition of H;
      if rank(H)==length(theta), go to 2
2. check rank of J
      if rank(J)<length(theta) - find out which parameters are involved
      if rank(J)==length(theta) => the parameters are identified at theta
by the moments included in J
