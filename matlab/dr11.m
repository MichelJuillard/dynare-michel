% Copyright (C) 2001 Michel Juillard
%
function dr=dr11(iorder,dr,cheik)

global M_ options_ oo_
global it_ stdexo_ means_ dr1_test_ bayestopt_

% hack for Bayes
global dr1_test_ bayestopt_

options_ = set_default_option(options_,'loglinear',0);

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv  = transpose(M_.lead_lag_incidence);
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;


if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('Error in model specification: some variables don"t appear as current') ;
end

if ~cheik
%  if xlen > 1
%    error (['SS: stochastic exogenous variables must appear only at the' ...
%      ' current period. Use additional endogenous variables']) ;
%  end
end
  
if M_.maximum_lead > 1 & iorder > 1
  error (['Models with leads on more than one period can only be solved' ...
	  ' at order 1'])
end

dr=set_state_space(dr);
kstate = dr.kstate;
kad = dr.kad;
kae = dr.kae;
nstatic = dr.nstatic;
nfwrd = dr.nfwrd;
npred = dr.npred;
nboth = dr.nboth;           
order_var = dr.order_var;
nd = size(kstate,1);

sdyn = M_.endo_nbr - nstatic;


tempex = oo_.exo_simul;

it_ = M_.maximum_lag + 1;
z = repmat(dr.ys,1,klen);
z = z(iyr0) ;
%M_.jacobia=real(diffext('ff1_',[z; oo_.exo_steady_state])) ;
%M_.jacobia=real(jacob_a('ff1_',[z; oo_.exo_steady_state])) ;
[junk,M_.jacobia] = feval([M_.fname '_dynamic'],z,oo_.exo_simul);
oo_.exo_simul = tempex ;
tempex = [];

nz = size(z,1);
k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_lag+1),:);
b = M_.jacobia(:,M_.lead_lag_incidence(M_.maximum_lag+1,order_var));
a = b\M_.jacobia(:,nonzeros(k1')); 
if any(isinf(a(:)))
  dr1_test_(1) = 5;
  dr1_test_(2) = bayestopt_.penalty;
end
if M_.exo_nbr
  fu = b\M_.jacobia(:,nz+1:end);
end

if M_.maximum_lead == 0 & M_.maximum_lag == 1;  % backward model with one lag
  dr.ghx = -a;
  dr.ghu = -fu;
  return;
elseif M_.maximum_lead == 0 & M_.maximum_lag > 1 % backward model with lags on more than
			       % one period
  e = zeros(endo_nbr,nd);				
  k = find(kstate(:,2) <= M_.maximum_lag+1 & kstate(:,4));
  e(:,k) = -a(:,kstate(k,4)) ;
  dr.ghx = e;
  dr.ghu = -fu;
end

% buildind D and E
d = zeros(nd,nd) ;
e = d ;

k = find(kstate(:,2) >= M_.maximum_lag+2 & kstate(:,3));
d(1:sdyn,k) = a(nstatic+1:end,kstate(k,3)) ;
k1 = find(kstate(:,2) == M_.maximum_lag+2);
a1 = eye(sdyn);
e(1:sdyn,k1) =  -a1(:,kstate(k1,1)-nstatic);
k = find(kstate(:,2) <= M_.maximum_lag+1 & kstate(:,4));
e(1:sdyn,k) = -a(nstatic+1:end,kstate(k,4)) ;
k2 = find(kstate(:,2) == M_.maximum_lag+1);
k2 = k2(~ismember(kstate(k2,1),kstate(k1,1)));
d(1:sdyn,k2) = a1(:,kstate(k2,1)-nstatic);

if ~isempty(kad)
  for j = 1:size(kad,1)
    d(sdyn+j,kad(j)) = 1 ;
    e(sdyn+j,kae(j)) = 1 ;
  end
end
options_ = set_default_option(options_,'qz_criterium',1.000001);
  
if  ~exist('mjdgges')
  % using Chris Sim's routines
  use_qzdiv = 1;
  [ss,tt,qq,w] = qz(e,d);
  [tt,ss,qq,w] = qzdiv(options_.qz_criterium,tt,ss,qq,w);
  ss1=diag(ss);
  tt1=diag(tt);
  warning_state = warning;
  warning off;
  oo_.eigenvalues = ss1./tt1 ;
  warning warning_state;
  nba = nnz(abs(eigval) > options_.qz_criterium);
else
  use_qzdiv = 0;
  [ss,tt,w,sdim,oo_.eigenvalues,info] = mjdgges(e,d,options_.qz_criterium);
  if info & info ~= nd+2;
    error(['ERROR' info ' in MJDGGES.DLL']);
  end
  nba = nd-sdim;
end

nyf = sum(kstate(:,2) > M_.maximum_lag+1);

if cheik
  dr.rank = rank(w(1:nyf,nd-nyf+1:end));
  % dr.eigval = oo_.eigenvalues;
  return
end

eigenvalues = sort(oo_.eigenvalues);

if nba > nyf;
%  disp('Instability !');  
  dr1_test_(1) = 3; %% More eigenvalues superior to unity than forward variables ==> instability.
  dr1_test_(2) = (abs(eigenvalues(nd-nba+1:nd-nyf))-1-1e-5)'*...
      (abs(eigenvalues(nd-nba+1:nd-nyf))-1-1e-5);% Distance to Blanchard-Khan conditions (penalty)
  return
elseif nba < nyf;
%  disp('Indeterminacy !');    
  dr1_test_(1) = 2; %% ==> Indeterminacy. 
  dr1_test_(2) = (abs(eigenvalues(nd-nyf+1:nd-nba))-1-1e-5)'*...
      (abs(eigenvalues(nd-nyf+1:nd-nba))-1-1e-5);% Distance to Blanchard-Khan conditions (penality)    
  %% warning('DR1: Blanchard-Kahn conditions are not satisfied. Run CHEIK to learn more!');
  return
end

np = nd - nyf;
n2 = np + 1;
n3 = nyf;
n4 = n3 + 1;
% derivatives with respect to dynamic state variables
% forward variables

if condest(w(1:n3,n2:nd)) > 1e9
%  disp('Indeterminacy !!');
  dr1_test_(1) = 2; 
  dr1_test_(2) = 1;
  return
end

warning_state = warning;
lastwarn('');
warning off;
gx = -w(1:n3,n2:nd)'\w(n4:nd,n2:nd)';

if length(lastwarn) > 0;
%  disp('Indeterminacy !!');
  dr1_test_(1) = 2; 
  dr1_test_(2) = 1;
  warning(warning_state);
  return
end

% predetermined variables
hx = w(1:n3,1:np)'*gx+w(n4:nd,1:np)';
hx = (tt(1:np,1:np)*hx)\(ss(1:np,1:np)*hx);

lastwarn('');
if length(lastwarn) > 0;
%  disp('Singularity problem in dr11.m');
  dr1_test_(1) = 2; 
  dr1_test_(2) = 1;
  warning(warning_state);
  return
end

k1 = find(kstate(n4:nd,2) == M_.maximum_lag+1);
k2 = find(kstate(1:n3,2) == M_.maximum_lag+2);
dr.ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];
  
%lead variables actually present in the model
j3 = nonzeros(kstate(:,3));
j4  = find(kstate(:,3));
% derivatives with respect to exogenous variables
if M_.exo_nbr
  a1 = eye(M_.endo_nbr);
  aa1 = [];
  if nstatic > 0
    aa1 = a1(:,1:nstatic);
  end
  dr.ghu = -[aa1 a(:,j3)*gx(j4,1:npred)+a1(:,nstatic+1:nstatic+ ...
						  npred) a1(:,nstatic+npred+1:end)]\fu;


    lastwarn('');
    if length(lastwarn) > 0;
%    disp('Singularity problem in dr11.m');
        dr1_test_(1) = 2; 
        dr1_test_(2) = 1;
        return
    end
end
warning(warning_state);

% static variables
if nstatic > 0
  temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
  j5 = find(kstate(n4:nd,4));
  temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
  dr.ghx = [temp; dr.ghx];
  temp = [];
end

if options_.loglinear == 1
    k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;

    dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
	     repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
    dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
end

% necessary when using Sims' routines
if use_qzdiv
  gx = real(gx);
  hx = real(hx);
  dr.ghx = real(dr.ghx);
  dr.ghu = real(dr.ghu);
end


