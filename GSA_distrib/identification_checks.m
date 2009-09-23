function [indok, indno, indweak] = identification_checks(H,JJ)
%In the  as2001 model, after H and J are computed,

% My suggestion is to have the following steps for identification check in
% dynare:

% 1. check rank of H at theta
npar = size(H,2);
if rank(H)<npar
  %         - find out which parameters are involved,
  % using something like the vnorm and the eigenvalue decomposition of H;
  ind1 = find(vnorm(H)~=0);
  if length(ind1)<npar,
    indno = find(~ismember([1:npar],ind1));
  end
  [e1,e2] = eig(H'*H);
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
  % near-dependencies, mostly due to strong pairwise colliniearities, which can
  % be checked using

  for ii = 1:size(H1,2);
    for jj = 1:size(H1,2);
      PcoH(ii,jj) = [cosn([H1(:,ii),H1(:,jj)])];
      PcoJ(ii,jj) = [cosn([JJ1(:,ii),JJ1(:,jj)])];
    end
  end
else % rank(H)==length(theta), go to 2
  % 2. check rank of J
  if rank(JJ)<size(JJ,2)
    %         - find out which parameters are involved
    ind1 = find(vnorm(JJ)~=0);
    if length(ind1)<npar,
      indno = find(~ismember([1:npar],ind1));
    end
    [e1,e2] = eig(J'*J);

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
%     near-dependencies, mostly due to strong pairwise colliniearities, which can
%     be checked using

    for ii = 1:size(H1,2);
      for jj = 1:size(H1,2);
        PcoH(ii,jj) = [cosn([H1(:,ii),H1(:,jj)])];
        PcoJ(ii,jj) = [cosn([JJ1(:,ii),JJ1(:,jj)])];
      end
    end
  else  %rank(J)==length(theta) =>
    disp('The parameters are identified at theta by the moments included in J')
  end
end






