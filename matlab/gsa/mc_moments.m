function [vdec, cc, ac] = mc_moments(mm, ss, dr)
global options_ M_

  [nr1, nc1, nsam] = size(mm);
  disp('Computing theoretical moments ...')
  h = waitbar(0,'Theoretical moments ...');
  
  for j=1:nsam,
    dr.ghx = mm(:, [1:(nc1-M_.exo_nbr)],j);
    dr.ghu = mm(:, [(nc1-M_.exo_nbr+1):end], j);
    if ~isempty(ss),
      set_shocks_param(ss(j,:));
    end
    [vdec(:,:,j), corr, autocorr, z, zz] = th_moments(dr,options_.varobs);
    cc(:,:,j)=triu(corr);
    dum=[];
    for i=1:options_.ar
    dum=[dum, autocorr{i}];
    end
    ac(:,:,j)=dum;
    waitbar(j/nsam,h)
  end
  close(h)
  disp(' ')
  disp('... done !')
