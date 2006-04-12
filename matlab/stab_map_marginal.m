function [proba, dproba] = stab_map_marginal(lpmat, ibehaviour, inonbehaviour, aname, ishock)
%function stab_map_1(lpmat, ibehaviour, inonbehaviour, aname, ishock)
%
% lpmat =  Monte Carlo matrix
% ibehaviour = index of behavioural runs
% inonbehaviour = index of non-behavioural runs
% ishock = 1 estimated shocks included
% ishock = 0 estimated shocks excluded (default)
%
% Plots: dotted lines for BEHAVIOURAL
%        solid lines for NON BEHAVIOURAL
% USES smirnov

global estim_params_ bayestopt_ M_ options_

if nargin<5,
    ishock=0;
end
fname_ = M_.fname;

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
%kernel_function = 'uniform';     % Gaussian kernel for Fast Fourrier Transform approximaton.  

if ishock,
    npar = nshock + estim_params_.np;
else
    npar = estim_params_.np;
end

for i=1:ceil(npar/12),
    figure,
    for j=1+12*(i-1):min(npar,12*i),
        subplot(3,4,j-12*(i-1))
        optimal_bandwidth = mh_optimal_bandwidth(lpmat(ibehaviour,j),length(ibehaviour),bandwidth,kernel_function); 
        [x1,f1] = kernel_density_estimate(lpmat(ibehaviour,j),number_of_grid_points,...
            optimal_bandwidth,kernel_function);
        plot(x1, f1,':k','linewidth',2)
        optimal_bandwidth = mh_optimal_bandwidth(lpmat(inonbehaviour,j),length(inonbehaviour),bandwidth,kernel_function); 
        [x1,f1] = kernel_density_estimate(lpmat(inonbehaviour,j),number_of_grid_points,...
            optimal_bandwidth,kernel_function);
        hold on, plot(x1, f1,'k','linewidth',2)
        
        %hist(lpmat(ibehaviour,j),30)
        if ishock,
            title(bayestopt_.name{j},'interpreter','none')
        else
            title(bayestopt_.name{j+nshock},'interpreter','none')
        end
    end
    saveas(gcf,[fname_,'_',aname,'_',int2str(i)])
    eval(['print -depsc2 ' fname_ '_' aname '_' int2str(i)]);
    eval(['print -dpdf ' fname_ '_' aname '_' int2str(i)]);
    if options_.nograph, close(gcf), end
end
