% Copyright (C) 2001 Michel Juillard
%
function dr=set_state_space(dr)

global M_ oo_ options_ it_

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('Error in model specification: some variables don"t appear as current') ;
end

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
nboth = length(both_var);
npred = length(pred_var);
nfwrd = length(fwrd_var);
nstatic = length(stat_var);
order_var = [ stat_var; pred_var; both_var; fwrd_var];
inv_order_var(order_var) = (1:M_.endo_nbr);

% building kmask for z state vector in t+1
if M_.maximum_lag > 0
  kmask = [];
  if M_.maximum_lead > 0 
    kmask = [cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+2:end,order_var)),1)] ;
  end
  kmask = [kmask; flipud(cumsum(M_.lead_lag_incidence(1:M_.maximum_lag,order_var),1))] ;
else
  kmask = cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+2:klen,order_var)),1) ;
end

kmask = kmask';
kmask = kmask(:);
i_kmask = find(kmask);          % index of nonzero entries in kmask
nd = size(i_kmask,1);           % size of the state vector
kmask(i_kmask) = [1:nd];

% auxiliary equations

% elements that are both in z(t+1) and z(t)
k1 = find([kmask(1:end-M_.endo_nbr) & kmask(M_.endo_nbr+1:end)] );
kad = [];
kae = [];
if ~isempty(k1)
  kad = kmask(k1+M_.endo_nbr);
  kae = kmask(k1);
end

% composition of state vector
% col 1: variable;           col 2: lead/lag in z(t+1); 
% col 3: A cols for t+1 (D); col 4: A cols for t (E)
kstate = [ repmat([1:M_.endo_nbr]',klen-1,1) kron([klen:-1:2]',ones(M_.endo_nbr,1)) ...
	   zeros((klen-1)*M_.endo_nbr,2)];
kiy = flipud(M_.lead_lag_incidence(:,order_var))';
kiy = kiy(:);
kstate(1:M_.maximum_lead*M_.endo_nbr,3) = kiy(1:M_.maximum_lead*M_.endo_nbr)-M_.endo_nbr;  
kstate(find(kstate(:,3) < 0),3) = 0;
kstate(M_.maximum_lead*M_.endo_nbr+1:end,4) = kiy((M_.maximum_lead+1)*M_.endo_nbr+1:end);  
% put in E only the current variables that are not already in D
kstate = kstate(i_kmask,:);

dr.order_var = order_var;
dr.inv_order_var = inv_order_var';
dr.nstatic = nstatic;
dr.npred = npred+nboth;
dr.kstate = kstate;
dr.kad = kad;
dr.kae = kae;
dr.nboth = nboth;
dr.nfwrd = nfwrd;
% number of forward variables in the state vector
dr.nsfwrd = sum(kstate(:,2) > M_.maximum_lag+1);
% number of predetermined variables in the state vector
dr.nspred = sum(kstate(:,2) <= M_.maximum_lag+1);

% copmutes column position of auxiliary variables for 
% compact transition matrix (only state variables)
aux = zeros(0,1);
k0 = kstate(find(kstate(:,2) <= M_.maximum_lag+1),:);;
i0 = find(k0(:,2) == M_.maximum_lag+1);
for i=M_.maximum_lag:-1:2
  i1 = find(k0(:,2) == i);
  n1 = size(i1,1);
  j = zeros(n1,1);
  for j1 = 1:n1
    j(j1) = find(k0(i0,1)==k0(i1(j1),1));
  end
  aux = [aux; i0(j)];
  i0 = i1;
end
dr.transition_auxiliary_variables = [(1:size(aux,1))' aux];
