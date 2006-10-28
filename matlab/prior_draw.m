function pdraw = prior_draw(init,cc)
% Build one draw from the prior distribution. 
% 
% INPUTS
%   o init           [integer]  scalar equal to 1 (first call) or 0.
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
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
global M_ options_ estim_params_  
persistent fname npar bounds pshape pmean pstd a b p3 p4 condition
  
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
        bounds = [-Inf Inf];
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
            b(i) = a*(1/mu - 1);        
          case 2;%Gamma prior
            mu = p1(i)-p3(i);
            b(i) = p2(i)^2/mu;
            a(i) = mu/b;
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
        pdraw(i) = rand*(p4(i)-p3(i)) + p3(i);
      case 3% Gaussian prior.
        while condition
            tmp = randn*pstd(i) + pmean(i);
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end
      case 2% Gamma prior.
        while condition
            g = gamma_draw(a(i),b(i),p3(i));
            if g >= bounds(i,1) && g <= bounds(i,2)
                pdraw(i) = g;
                break
            end
        end
      case 1% Beta distribution (TODO: generalized beta distribution)
        while condition
            y1 = gamma_draw(a(i),1,0);
            y2 = gamma_draw(b(i),1,0);
            tmp = y1/(y1+y2);
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = pmean(i)+tmp*pstd(i);
                break
            end
        end
      case 4% INV-GAMMA1 distribution
        while condition
            tmp = sqrt(1/gamma_draw(p2(i)/2,1/p1(i),0));
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end        
      case 6% INV-GAMMA2 distribution  
        while condition
            tmp = 1/gamma_draw(p2(i)/2,1/p1(i),0);
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end
      otherwise
        % Nothing to do here.
    end
end

  

  
  
function gamma_draw(a,b,c)
% Bauwens, Lubrano & Richard (page 316)
if a >30
    z = randn;
    g = b*(z+sqrt(4*a-1))^2/4 + c; 
else
    condi = 1
    while condi
        x = -1;
        while x<0
            u1 = rand;
            y = tan(pi*u1);
            x = y*sqrt(2*a-1)+a-1; 
        end
        u2 = rand;
        if log(u2) <= log(1+y^2)+(a-1)*log(x/(a-1))-y*sqrt(2*a-1);
            break
        end
    end
    g = x*b+c;
end