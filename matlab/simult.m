% Copyright (C) 2001 Michel Juillard
%

function y_=simult(ys, dr)
global M_ options_ oo_
global  it_ means_

order = options_.order;
replic = options_.replic;
if replic == 0
  replic = 1;
end
seed = options_.simul_seed;
options_.periods = options_.periods;

it_ = M_.maximum_lag + 1 ;

if replic > 1
  fname = [M_.fname,'_simul'];
  fh = fopen(fname,'w+');
end

% eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0));
nxs = length(i_exo_var);
oo_.exo_simul = zeros(M_.maximum_lag+M_.maximum_lead+options_.periods,M_.exo_nbr);
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));

for i=1:replic
  if isempty(seed)
    randn('state',sum(100*clock));
  else
    randn('state',seed+i-1);
  end
  if ~isempty(M_.Sigma_e)
    oo_.exo_simul(:,i_exo_var) = randn(M_.maximum_lag+M_.maximum_lead+options_.periods,nxs)*chol_S;
  end
  y_ = simult_(ys,dr,oo_.exo_simul,order);
  if replic > 1
    fwrite(fh,oo_.endo_simul(:,M_.maximum_lag+1:end),'float64');
  end
end

if replic > 1
  fclose(fh);
end


% 02/20/01 MJ replaced ys by dr.ys
% 02/22/01 MJ removed commented out lines
%             removed useless temps
%             stderr_ replaced by M_.Sigma_e
% 02/28/01 MJ changed expression for M_.Sigma_e
% 02/18/03 MJ added ys in the calling sequence for arbitrary initial values
%             suppressed useless calling parameter istoch
% 05/10/03 MJ removed repmat() in call to simult_() for lag > 1
% 05/29/03 MJ test for 0 variances
