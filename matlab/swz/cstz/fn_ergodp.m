function gpi = fn_ergodp(P)
% gpi = fn_ergodp(P)
%    Compute the ergodic probabilities.  See Hamilton p.681.
%
% P:  n-by-n matrix of transition matrix where all elements in each column sum up to 1.
%-----
% gpi:  n-by-1 vector of ergodic probabilities.
%
% Tao Zha August 2000

[gpim,gpid] = eig(P);  % m: matrix; d: diagonal
[gpidv,gpidvinx] = sort(diag(gpid));
gpidv = fliplr(gpidv);
gpidvinx = flipud(gpidvinx);
gpim = gpim(:,gpidvinx);
gpi = gpim(:,1)/sum(gpim(:,1));
