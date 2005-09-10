function [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(xx,info)
% stephane.adjemian@ens.fr [09-09-2005]
% Computes posterior mean, median, variance, HPD interval, deciles, and density from posterior draws.
global options_

xx = xx(:);
xx = sort(xx);


post_mean = mean(xx);
post_median = median(xx);
post_var = var(xx);

n = length(xx);
m = round((1-options_.mh_conf_sig)*n);
k = zeros(m,1);
jj = n-m;
for ii = 1:m
  k(ii) = xx(jj)-xx(ii);
  jj = jj + 1;
end
[kmin,idx] = min(k);
hpd_interval = [xx(idx) xx(idx)+kmin];

post_deciles = xx([round(0.1*n) ...
       round(0.2*n)...
       round(0.3*n)...
       round(0.4*n)...
       round(0.5*n)...
       round(0.6*n)...
       round(0.7*n)...
       round(0.8*n)...
       round(0.9*n)]);

if ~info
  density = [];
else
  number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
  bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
  kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
  optimal_bandwidth = mh_optimal_bandwidth(xx,length(xx),bandwidth,kernel_function);
  % [abscissa,f] = kernel_density_estimate(,,,)
  [density(:,1),density(:,2)] = kernel_density_estimate(xx,number_of_grid_points,optimal_bandwidth,kernel_function);
end