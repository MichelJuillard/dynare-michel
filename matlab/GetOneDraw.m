function [xparams, logpost] = GetOneDraw(type)
% stephane.adjemian@ens.fr [09-25-2005]

  switch type
   case 'posterior'
    [xparams, logpost] = metropolis_draw(0);
   case 'prior'
    xparams = prior_draw(0);
  end