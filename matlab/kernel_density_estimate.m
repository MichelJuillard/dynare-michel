function [abscissa,f] = kernel_density_estimate(data,number_of_grid_points,bandwidth,kernel_function) 

% function [abscissa,f] = kernel_density_estimate(data,number_of_grid_points,bandwidth,kernel_function) 
% Estimates a continuous density. 
% 
% INPUTS
%    data:                    data
%    number_of_grid_points:   number of grid points
%    bandwidth:               scalar equals to 0,-1 or -2. For a value different from 0,-1 or -2 the
%                             function will return optimal_bandwidth = bandwidth.
%    kernel_function:         'gaussian','uniform','triangle','epanechnikov',
%                             'quartic','triweight','cosinus'.
%
% OUTPUTS
%    abscissa:                value on the abscissa axis
%    f:                       density
%        
% SPECIAL REQUIREMENTS
%    A kernel density estimator is used (see Silverman [1986], "Density estimation for statistics and data analysis")
%    The code is adapted from Anders Holtsberg's matlab toolbox (stixbox). 
%
% part of DYNARE, copyright Dynare Team (2004-2008)
% Gnu Public License.



if size(data,2) > 1 & size(data,1) == 1 
    data = transpose(data); 
elseif size(data,2)>1 & size(data,1)>1 
    error('kernel_density_estimate :: data must be a one dimensional array.'); 
end
test = log(number_of_grid_points)/log(2);
if (abs(test-round(test)) > 10^(-12))
    error('kernel_density_estimate :: The number of grid points must be a power of 2.');
end

n = size(data,1);

%% KERNEL SPECIFICATION...
if strcmpi(kernel_function,'gaussian') 
    k    = inline('inv(sqrt(2*pi))*exp(-0.5*x.^2)'); 
elseif strcmpi(kernel_function,'uniform') 
    k    = inline('0.5*(abs(x) <= 1)'); 
elseif strcmpi(kernel_function,'triangle') 
    k    = inline('(1-abs(x)).*(abs(x) <= 1)'); 
elseif strcmpi(kernel_function,'epanechnikov') 
    k    = inline('0.75*(1-x.^2).*(abs(x) <= 1)'); 
elseif strcmpi(kernel_function,'quartic') 
    k    = inline('0.9375*((1-x.^2).^2).*(abs(x) <= 1)'); 
elseif strcmpi(kernel_function,'triweight') 
    k    = inline('1.09375*((1-x.^2).^3).*(abs(x) <= 1)'); 
elseif strcmpi(kernel_function,'cosinus') 
    k    = inline('(pi/4)*cos((pi/2)*x).*(abs(x) <= 1)'); 
end

%% COMPUTE DENSITY ESTIMATE... Gaussian kernel should be used (FFT).
a  = min(data) - (max(data)-min(data))/3;
b  = max(data) + (max(data)-min(data))/3;
abscissa = linspace(a,b,number_of_grid_points)';
d  = abscissa(2)-abscissa(1); 
xi = zeros(number_of_grid_points,1);
xa = (data-a)/(b-a)*number_of_grid_points; 
for i = 1:n;
    indx = floor(xa(i));
    temp = xa(i)-indx;
    xi(indx+[1 2]) = xi(indx+[1 2]) + [1-temp,temp]';
end;
xk = [-number_of_grid_points:number_of_grid_points-1]'*d;
kk = k(xk/bandwidth);
kk = kk / (sum(kk)*d*n);
f = ifft(fft(fftshift(kk)).*fft([xi ;zeros(size(xi))]));
f = real(f(1:number_of_grid_points));