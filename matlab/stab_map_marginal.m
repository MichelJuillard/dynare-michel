function stab_map_marginal(lpmat, ibehaviour, inonbehaviour, aname, ipar, dirname)
%function stab_map_1(lpmat, ibehaviour, inonbehaviour, aname, ipar, dirname)
%
% lpmat =  Monte Carlo matrix
% ibehaviour = index of behavioural runs
% inonbehaviour = index of non-behavioural runs
% aname = label of the analysis
% ipar = index array of parameters to plot
% dirname = (OPTIONAL) path of the directory where to save 
%            (default: current directory)
%
% Plots: dotted lines for BEHAVIOURAL
%        solid lines for NON BEHAVIOURAL
% USES smirnov

global estim_params_ bayestopt_ M_ options_

fname_ = M_.fname;
if nargin<5,
  ipar=[1:npar];
end
nparplot=length(ipar);
if nargin<6,
  dirname='';;
end

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

npar=size(lpmat,2);
ishock= npar>estim_params_.np;

number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
%kernel_function = 'uniform';     % Gaussian kernel for Fast Fourrier Transform approximaton.  

lpmat=lpmat(:,ipar);
ftit=bayestopt_.name(ipar+nshock*(1-ishock));

for i=1:ceil(npar/12),
    figure('name',aname),
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
        
        title(ftit{j},'interpreter','none')
    end
    saveas(gcf,[dirname,'/',fname_,'_',aname,'_',int2str(i)])
    eval(['print -depsc2 ' dirname '\' fname_ '_' aname '_' int2str(i)]);
    eval(['print -dpdf ' dirname '\' fname_ '_' aname '_' int2str(i)]);
    if options_.nograph, close(gcf), end
end
