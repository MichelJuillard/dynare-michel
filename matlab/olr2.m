% Copyright (C) 2003 Michel Juillard
%
% computes an optimal policy as the optimal linear regulator
%
function dr=olr2(dr,olr_inst,bet,obj_var,W)

global M_ options_ oo_
global it_ stdexo_

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;

inst_nbr = size(olr_inst,1);
inst_i = zeros(inst_nbr,1);
for i=1:inst_nbr
  k = strmatch(olr_inst(i,:),M_.endo_names,'exact');
  if isempty(k)
    error(sprintf('OLR_INST %s isn''t a declared variable'));
  else
    inst_i(i) = k;
  end
end

if M_.maximum_lead == 0
  error ('OLR : No forward variable: no point in using OLR') ;
end

if find(any(M_.lead_lag_incidence([1:M_.maximum_lag M_.maximum_lag+2:M_.maximum_lead],inst_i),2))
  error('OLR: instruments can only appear at the current period');
end

non_inst_i = setdiff([1:M_.endo_nbr],inst_i);
iy1_ = M_.lead_lag_incidence(:,non_inst_i);
endo_nbr_1 = M_.endo_nbr - inst_nbr;

if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('Error in model specification: some variables don''t appear as current') ;
end

if xlen > 1
  error (['SS: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end

tempex = oo_.exo_simul;

it_ = M_.maximum_lag + 1;
dr.ys = oo_.steady_state; 
z = repmat(dr.ys,1,klen);
z = z(iyr0) ;

%M_.jacobia=real(jacob_a('ff1_',[z; oo_.exo_steady_state])) ;
[junk,M_.jacobia] = feval([M_.fname '_dynamic'],z,oo_.exo_simul);

oo_.exo_simul = tempex ;
tempex = [];

nz = size(z,1);
fwrd_var = find(any(M_.lead_lag_incidence(M_.maximum_lag+2:end,:),1))';
if M_.maximum_lag > 0
  pred_var = find(any(M_.lead_lag_incidence(1:M_.maximum_lag,:),1))';
  both_var = intersect(pred_var,fwrd_var);
  pred_var = setdiff(pred_var,both_var);
  fwrd_var = setdiff(fwrd_var,both_var);
  stat_var = setdiff([1:M_.endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
  pred_var = [];
  both_var = [];
  stat_var = setdiff([1:M_.endo_nbr]',fwrd_var);
end

stat_var = setdiff(stat_var,inst_i);
% static variables in objective function
[stat_obj_var] = intersect(obj_var,stat_var);
n_stat_obj_var = length(stat_obj_var);
pred_var = [stat_obj_var; pred_var];
nboth = length(both_var);
npred = length(pred_var);
nfwrd = length(fwrd_var);
nstatic1 = length(stat_var);
nstatic = nstatic1-n_stat_obj_var;

order_var = [ setdiff(stat_var,obj_var); pred_var; both_var; fwrd_var];
k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_lag+1),:);
b = M_.jacobia(1:endo_nbr_1,M_.lead_lag_incidence(M_.maximum_lag+1,order_var));
a = b\M_.jacobia(1:endo_nbr_1,nonzeros(k1')); 
if M_.exo_nbr
  fu = -b\M_.jacobia(1:endo_nbr_1,nz+1:end);
end
% instruments' effects
b = -b\M_.jacobia(1:endo_nbr_1,M_.lead_lag_incidence(M_.maximum_lag+1,inst_i));
% building kmask for z state vector in t+1
if M_.maximum_lag > 0
  if M_.maximum_lead > 0 
    kmask = [flipud(cumsum(M_.lead_lag_incidence(1:M_.maximum_lag,order_var),1))] ;
    kmask = [kmask; flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+2:end,order_var)),1))] ;
  end
else
  kmask = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+2:klen,order_var)),1));
end
% static variables in objective function
kmask(1,nstatic+1:nstatic1) = 1;
kmask = kmask';

% lags aren't ordered as in dr1 !
% this is necessary to have the zeros on R's diagonal aligned with jump variables
kmask = kmask(:);
i_kmask = find(kmask);          % index of nonzero entries in kmask
nd = size(i_kmask,1);           % size of the state vector
kmask(i_kmask) = [1:nd];

% auxiliary equations

% elements that are both in z(t+1) and z(t)
kad = [];
kae = [];
k1 = find([kmask(1:endo_nbr_1) & kmask(M_.maximum_lag*endo_nbr_1+1:(M_.maximum_lag+1)*endo_nbr_1)] );
if ~isempty(k1)
  kad = kmask(k1);
  kae = kmask(k1+M_.maximum_lag*endo_nbr_1);
end

if M_.maximum_lag > 1
  k1 = find([kmask(endo_nbr_1+1:M_.maximum_lag*endo_nbr_1) & kmask(1:(M_.maximum_lag-1)*endo_nbr_1)] );
  if ~isempty(k1)
    kad = [kad; kmask(k1+endo_nbr_1)];
    kae = [kae; kmask(k1)];
  end
end
if M_.maximum_lead > 1
  k1 = find([kmask((M_.maximum_lag+1)*endo_nbr_1+1:end) & kmask(M_.maximum_lag*endo_nbr_1+1:end-endo_nbr_1)] );
  if ~isempty(k1)
    kad = [kad; kmask(M_.maximum_lag*endo_nbr_1+k1)];
    kae = [kae; kmask((M_.maximum_lag+1)*endo_nbr_1+k1)];
  end
end

% composition of state vector
% col 1: variable;           col 2: lead/lag in z(t+1); 
% col 3: A cols for t+1 (D); col 4: A cols for t (E)
kstate = [ repmat([1:endo_nbr_1]',klen-1,1) ...
	   [kron([M_.maximum_lag+1:-1:2]',ones(endo_nbr_1,1)); ...
	    kron([M_.maximum_lag+2:klen]',ones(endo_nbr_1,1))] ...
	   zeros((klen-1)*endo_nbr_1,2)];
kiy = [flipud(M_.lead_lag_incidence(1:M_.maximum_lag+1,order_var)); M_.lead_lag_incidence(M_.maximum_lag+2:end,order_var)]';
kiy = kiy(:);
i1 = find(kiy((M_.maximum_lag+1)*endo_nbr_1+1:end));
kstate(i1+M_.maximum_lag*endo_nbr_1,3) = kiy(i1+(M_.maximum_lag+1)*endo_nbr_1)-M_.endo_nbr;  
kstate(1:M_.maximum_lag*endo_nbr_1,4) = kiy(endo_nbr_1+1:(M_.maximum_lag+1)*endo_nbr_1);  
% put in E only the current variables that are not already in D
kstate = kstate(i_kmask,:);

sdyn = M_.endo_nbr-nstatic-inst_nbr;

% buildind D and E
d = zeros(nd,nd) ;
e = d ;

% dynamical system
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
[Q,R] = qr(d);
C = Q'*[e [[b(nstatic+1:end,:) fu(nstatic+1:end,:)];...
	   zeros(nd-sdyn,inst_nbr+M_.exo_nbr)]];

nyf = sum(kstate(:,2) > M_.maximum_lag+1);


np = nd - nyf;


dd = zeros(2*nd+inst_nbr,2*nd+inst_nbr);
ee = dd;

dd(1:nd,1:np) = R(:,1:np);
ee(1:nd,1:np) = C(:,1:np);
dd(1:nd,nd+1:nd+nyf) = R(:,np+1:end);
ee(1:nd,nd+1:nd+nyf+inst_nbr) = C(:,np+1:nd+inst_nbr);


% derivatives with respect to x
% building QQ and UU
m = 1;
v0 = kstate(find(kstate(:,2)==M_.maximum_lag+2),1);
for i=1:nd;
  if (kstate(i,2) == M_.maximum_lag+1 & isempty(find(kstate(i,1)==v0)))
    k = find(order_var(kstate(i,1))==obj_var);
    if ~isempty(k)
      iq(m) = i;
      m = m+1;
    end
  elseif kstate(i,2) == M_.maximum_lag+2
    k = find(order_var(kstate(i,1))==obj_var);
    if ~isempty(k)
      iq(m) = i;
      m = m+1;
    end
  end
end
QQ1 = zeros(nd,nd);
QQ1(iq,iq) = W(obj_var,obj_var);
UU1 = zeros(nd,inst_nbr);
UU1(iq,:) = W(obj_var,inst_i);
RR = W(inst_i,inst_i);
offset = nd;
ee(offset+1:2*offset,1:np) = bet*QQ1(1:np,:)';
dd(offset+1:2*offset,np+1:nd) = bet*C(np+1:nd,1:nd)';
ee(offset+1:2*offset,nd+1:nd+nyf) = bet*QQ1(np+1:end,:)';
dd(offset+1:2*offset,nd+nyf+inst_nbr+1:end)=bet*C(1:np,1:nd)';
ee(offset+1:2*offset,np+1:nd) = R(np+1:end,:)';
ee(offset+1:2*offset,nd+nyf+1:nd+nyf+inst_nbr) = -bet*UU1;
ee(offset+1:2*offset,nd+nyf+inst_nbr+1:end) = R(1:np,:)';

%derivatives with respect to u
offset = 2*nd;
dd(offset+1:end,np+1:nd) = -C(np+1:end,nd+1:nd+inst_nbr)';
dd(offset+1:end,nd+nyf+inst_nbr+1:end) = -C(1:np,nd+1:nd+inst_nbr)';
ee(offset+1:end,1:np) = UU1(1:np,:)';
ee(offset+1:end,nd+1:nd+nyf) = UU1(np+1:end,:)';
ee(offset+1:end,nd+nyf+1:nd+nyf+inst_nbr) = RR;

[ss,tt,w,sdim,oo_.eigenvalues,info] = mjdgges(ee,dd);

nba = sum(abs(oo_.eigenvalues) > 1+1e-5);
nyf1 = nd+inst_nbr;
nd = 2*nd+inst_nbr;

if nba ~= nyf1;
  warning('OLR: Blanchard-Kahn conditions are not satisfied.');
end

np1 = nd - nyf1;
n2 = np1 + 1;
n3 = np1;
n4 = n3 + 1;

% derivatives with respect to dynamic state variables
% forward variables
gx = -w(n4:nd,n2:nd)'\w(1:n3,n2:nd)';
% predetermined variables
hx = w(n4:nd,1:np1)'*gx+w(1:n3,1:np1)';
hx = (tt(1:np1,1:np1)*hx)\(ss(1:np1,1:np1)*hx);

% including Lagrange multipliers in M_.endo_names, order_var and kstate
for i=1:M_.maximum_lead
  k = find(kstate(:,2)==M_.maximum_lag+1+i);
  temp = strcat('mult_',M_.endo_names(order_var(kstate(k,1)),:));
  temp = strcat(temp,['_' int2str(i)]);
  M_.endo_names = strvcat(M_.endo_names,temp);
end

if nyf > nboth	
  order_var = [order_var(1:nstatic+npred+nboth);[M_.endo_nbr+1:M_.endo_nbr+nyf]'; ...
		  order_var(nstatic+npred+nboth+1:end); inst_i];
else
  order_var = [order_var(1:nstatic+npred+nboth);[M_.endo_nbr+1:M_.endo_nbr+nyf]'; ...
		  inst_i];
end
kstate = [kstate(1:np,:);zeros(nyf,4);kstate(np+1:end,:);zeros(inst_nbr+np,4)];
k = find(kstate(:,2) <= M_.maximum_lag+1);
kstate(np+1:np+nyf,1:2) = [[nstatic+npred+nboth+1:nstatic+npred+nboth+nyf]' ...
		    (M_.maximum_lag+1)*ones(nyf,1)];
kstate(np+2*nyf+1:np+2*nyf+inst_nbr,1:2) = [[endo_nbr_1+1:endo_nbr_1+inst_nbr]' ...
		    (M_.maximum_lag+2)*ones(inst_nbr,1)];
kstate(np+2*nyf+inst_nbr+1:end,1:2) = [[M_.endo_nbr+1:M_.endo_nbr+np]' (M_.maximum_lag+2)*ones(np,1)];

k1 = find(kstate(:,2) == M_.maximum_lag+1);
k2 = find(kstate(:,2) == M_.maximum_lag+2)-np-nyf;
dr.ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];

%lead variables actually present in the model
% derivatives with respect to exogenous variables
if M_.exo_nbr
  n1 = find(kstate(:,2) > M_.maximum_lag+1);
  ghu = [dd(:,n1(1):end)*gx+dd(:,1:n1(1)-1) -ee(:,np+nyf+1:end)]\[C(:,end-M_.exo_nbr+1:end); zeros(size(dd,1)-size(C,1),M_.exo_nbr)];
  dr.ghu = [ghu(k1,:);ghu(k2(nboth+1:end)+np+nyf,:)];
end

% static variables
if nstatic > 0
  j3 = nonzeros(kstate(:,3));
  j4  = find(kstate(:,3))-np-nyf;
  temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
  temp = temp + b(1:nstatic,:)*gx(nyf+1:nyf+inst_nbr,:);
  j5 = find(kstate(:,4));
  temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
  dr.ghx = [temp; dr.ghx];
  temp = -a(1:nstatic,j3)*gx(j4,:)*ghu(1:np+nyf,:);
  temp = temp + b(1:nstatic,:)*ghu(np+2*nyf+1:np+2*nyf+inst_nbr,:);
  temp = temp + fu(1:nstatic,:);
  dr.ghu = [temp; dr.ghu]; 
  temp = [];
end

dr.ghx = dr.ghx(1:M_.endo_nbr+nyf,:);
dr.ghu = dr.ghu(1:M_.endo_nbr+nyf,:);

dr.ys = [dr.ys; zeros(nyf,1)];
dr.nstatic1 = nstatic1;
dr.nstatic = nstatic;
dr.npred = npred+nyf+nboth;
dr.kstate = kstate;
dr.order_var = order_var;
M_.endo_nbr = M_.endo_nbr+nyf;

% 05/29/03 MJ replaced diffext by jacobia (much faster)
%             corrected kmask for static variables in objective function




