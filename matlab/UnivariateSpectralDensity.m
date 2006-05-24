% stephane.adjemian@ens.fr
function [omega,f] = UnivariateSpectralDensity(dr,var_list)
% This function computes the theoretical spectral density of each
% endogenous variable declared in var_list. Results are stored in 
% oo_ and may be plotted.
% 
% Adapted from th_autocovariances.m.  
global options_ oo_ M_

omega = []; f = [];

if options_.order > 1
  disp('UnivariateSpectralDensity :: I Cannot compute the theoretical spectral density') 
  disp('with a second order approximation of the DSGE model!')
  disp('Set order = 1.')
  return
end

pltinfo  = 1;%options_.SpectralDensity.Th.plot;
cutoff   = 100;%options_.SpectralDensity.Th.cutoff;
sdl      = 0.1;%options_.SepctralDensity.Th.sdl;
omega    = (0:sdl:pi)';
GridSize = length(omega);
exo_names_orig_ord  = M_.exo_names_orig_ord;
if sscanf(version('-release'),'%d') < 13
  warning off
else
  eval('warning off MATLAB:dividebyzero')
end
if nargin<2
  var_list = [];
end
nvar = size(var_list,1);
if nvar == 0
  nvar = length(dr.order_var);
  ivar = [1:nvar]';
else
    ivar=zeros(nvar,1);
    for i=1:nvar
      i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
      	error (['One of the variable specified does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
end
f = zeros(nvar,GridSize);
ghx = dr.ghx;
ghu = dr.ghu;
npred = dr.npred;
nstatic = dr.nstatic;
kstate = dr.kstate;
order = dr.order_var;
iv(order) = [1:length(order)];
nx = size(ghx,2);
ikx = [nstatic+1:nstatic+npred];
A = zeros(nx,nx);
k0 = kstate(find(kstate(:,2) <= M_.maximum_lag+1),:);
i0 = find(k0(:,2) == M_.maximum_lag+1);
i00 = i0;
n0 = length(i0);
A(i0,:) = ghx(ikx,:);
AS = ghx(:,i0);
ghu1 = zeros(nx,M_.exo_nbr);
ghu1(i0,:) = ghu(ikx,:);
for i=M_.maximum_lag:-1:2
  i1 = find(k0(:,2) == i);
  n1 = size(i1,1);
  j1 = zeros(n1,1);
  j2 = j1;
  for k1 = 1:n1
    j1(k1) = find(k0(i00,1)==k0(i1(k1),1));
    j2(k1) = find(k0(i0,1)==k0(i1(k1),1));
  end
  AS(:,j1) = AS(:,j1)+ghx(:,i1);
  i0 = i1;
end
b = ghu1*M_.Sigma_e*ghu1';
[A,B] = kalman_transition_matrix(dr);
% index of predetermined variables in A
i_pred = [nstatic+(1:npred) M_.endo_nbr+1:length(A)];
A = A(i_pred,i_pred);
[vx, ns_var] =  lyapunov_symm(A,b);
i_ivar = find(~ismember(ivar,dr.order_var(ns_var+nstatic)));
ivar = ivar(i_ivar);
iky = iv(ivar);  
aa = ghx(iky,:);
bb = ghu(iky,:);
Gamma = zeros(nvar,cutoff+1);
tmp = aa*vx*aa'+ bb*M_.Sigma_e*bb';
k = find(abs(tmp) < 1e-12);
tmp(k) = 0;
Gamma(:,1) = diag(tmp);
vxy = (A*vx*aa'+ghu1*M_.Sigma_e*bb');
tmp = aa*vxy;
k = find(abs(tmp) < 1e-12);
tmp(k) = 0;
Gamma(:,2) = diag(tmp);
for i=2:cutoff
  vxy = A*vxy;
  tmp = aa*vxy;
  k = find(abs(tmp) < 1e-12);
  tmp(k) = 0;
  Gamma(:,i+1) = diag(tmp);
end
H = 1:cutoff;
for i=1:nvar
  f(i,:) = Gamma(i,1)/(2*pi) + Gamma(i,H+1)*cos(H'*omega')/pi;
end  

if sscanf(version('-release'),'%d') < 13
  warning on
else
  eval('warning on MATLAB:dividebyzero')
end

if pltinfo
  for i= 1:nvar
    figure('Name',['Spectral Density of ' deblank(M_.endo_names(ivar(i),:)) '.'])
    plot(omega,f(i,:),'-k','linewidth',2)
    xlabel('0 \leq \omega \leq \pi')
    ylabel('f(\omega)')
    box on
    axis tight
  end
end