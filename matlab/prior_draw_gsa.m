function pdraw = prior_draw_gsa(init,rdraw,cc)
% Draws from the prior distributions 
% Adapted by M. Ratto from prior_draw (of DYNARE, copyright M. Juillard), 
% for use with Sensitivity Toolbox for DYNARE
% 
% 
% INPUTS
%   o init           [integer]  scalar equal to 1 (first call) or 0.
%   o rdraw          
%   o cc             [double]   two columns matrix (same as in
%                               metropolis.m), constraints over the
%                               parameter space (upper and lower bounds).
%    
% OUTPUTS 
%   o pdraw          [double]   draw from the joint prior density. 
%
% ALGORITHM 
%   ...       
%
% SPECIAL REQUIREMENTS
%   MATLAB Statistics Toolbox
%  
%  
% Part of the Sensitivity Analysis Toolbox for DYNARE
%
% Written by Marco Ratto, 2006
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Disclaimer: This software is not subject to copyright protection and is in the public domain. 
% It is an experimental system. The Joint Research Centre of European Commission 
% assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.
%

global M_ options_ estim_params_  bayestopt_
persistent fname npar bounds pshape pmean pstd a b p1 p2 p3 p4 condition
  
if init
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    MhDirectoryName = CheckPath('metropolis');
    fname = [ MhDirectoryName '/' M_.fname];
    pshape = bayestopt_.pshape;
    pmean = bayestopt_.pmean;
    pstd  = bayestopt_.pstdev;
    p1 = bayestopt_.p1;
    p2 = bayestopt_.p2;
    p3 = bayestopt_.p3;
    p4 = bayestopt_.p4;
    a = zeros(npar,1);
    b = zeros(npar,1); 
    if nargin == 2
        bounds = cc;
    else
        bounds = kron(ones(npar,1),[-Inf Inf]);
    end
    for i = 1:npar
        switch pshape(i)
          case 3% Gaussian prior
            b(i) = pstd(i)^2/(pmean(i)-p3(i));
            a(i) = (pmean(i)-p3(i))/b(i);
          case 1% Beta prior
            mu = (p1(i)-p3(i))/(p4(i)-p3(i));
            stdd = p2(i)/(p4(i)-p3(i));
            a(i) = (1-mu)*mu^2/stdd^2 - mu;
            b(i) = a(i)*(1/mu - 1);        
          case 2;%Gamma prior
            mu = p1(i)-p3(i);
            b(i) = p2(i)^2/mu;
            a(i) = mu/b(i);
          case {5,4,6}
            % Nothing to do here
            %
            % 4: Inverse gamma, type 1, prior
            %    p2(i) = nu
            %    p1(i) = s
            % 6: Inverse gamma, type 2, prior
            %    p2(i) = nu
            %    p1(i) = s
            % 5: Uniform prior
            %    p3(i) and p4(i) are used.
          otherwise
            disp('prior_draw :: Error!')
            disp('Unknown prior shape.')
            return
        end
        pdraw = zeros(npar,1);  
    end
    condition = 1;
    pdraw = zeros(npar,1);
    return
end


for i = 1:npar
    switch pshape(i)
      case 5% Uniform prior.
        pdraw(:,i) = rdraw(:,i)*(p4(i)-p3(i)) + p3(i);
      case 3% Gaussian prior.
        pdraw(:,i) = norminv(rdraw(:,i),pmean(i),pstd(i));
      case 2% Gamma prior.
        pdraw(:,i) = gaminv(rdraw(:,i),a(i),b(i))+p3(i);
      case 1% Beta distribution (TODO: generalized beta distribution)
        pdraw(:,i) = betainv(rdraw(:,i),a(i),b(i))*(p4(i)-p3(i))+p3(i);
      case 4% INV-GAMMA1 distribution 
        % TO BE CHECKED
        pdraw(:,i) =  sqrt(1./gaminv(rdraw(:,i),p2(i)/2,1/p1(i)));
      case 6% INV-GAMMA2 distribution  
        % TO BE CHECKED
        pdraw(:,i) =  1./gaminv(rdraw(:,i),p2(i)/2,1/p1(i));
      otherwise
        % Nothing to do here.
    end
end

  
