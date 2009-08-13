function pdraw = prior_draw_gsa(init,rdraw)
% Draws from the prior distributions 
% Adapted by M. Ratto from prior_draw (of DYNARE, copyright M. Juillard), 
% for use with Sensitivity Toolbox for DYNARE
% 
% 
% INPUTS
%   o init           [integer]  scalar equal to 1 (first call) or 0.
%   o rdraw          
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

% global M_ options_ estim_params_  bayestopt_
global bayestopt_
persistent npar pshape p6 p7 p3 p4
  
if init
    pshape = bayestopt_.pshape;
    p6 = bayestopt_.p6;
    p7 = bayestopt_.p7;
    p3 = bayestopt_.p3;
    p4 = bayestopt_.p4;
    npar = length(p6);
    pdraw = zeros(npar,1);
    return
end


for i = 1:npar
    switch pshape(i)
      case 5% Uniform prior.
        pdraw(:,i) = rdraw(:,i)*(p4(i)-p3(i)) + p3(i);
      case 3% Gaussian prior.
        pdraw(:,i) = norm_inv(rdraw(:,i),p6(i),p7(i));
      case 2% Gamma prior.
        pdraw(:,i) = gamm_inv(rdraw(:,i),p6(i),p7(i))+p3(i);
      case 1% Beta distribution (TODO: generalized beta distribution)
        pdraw(:,i) = beta_inv(rdraw(:,i),p6(i),p7(i))*(p4(i)-p3(i))+p3(i);
      case 4% INV-GAMMA1 distribution 
        % TO BE CHECKED
        pdraw(:,i) =  sqrt(1./gamm_inv(rdraw(:,i),p7(i)/2,2/p6(i)))+p3(i);
      case 6% INV-GAMMA2 distribution  
        % TO BE CHECKED
        pdraw(:,i) =  1./gamm_inv(rdraw(:,i),p7(i)/2,2/p6(i))+p3(i);
      otherwise
        % Nothing to do here.
    end
end

  
