function err = ig1fun(nu,mu2,sigma2)
    err = 2*mu2*gamma(nu/2).^2-(sigma2+mu2)*(nu-2).*gamma((nu-1)/2).^2;