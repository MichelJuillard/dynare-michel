% Copyright (C) 2001 Michel Juillard
%
function [dr,info]=dr1(dr,task)
global M_ options_ oo_

global olr_state
% info = 1: the model doesn't define current variables uniquely
% info = 2: problem in mjdgges.dll info(2) contains error code
% info = 3: BK order condition not satisfied info(2) contains "distance"
%           absence of stable trajectory
% info = 4: BK order condition not satisfied info(2) contains "distance"
%           indeterminacy
% info = 5: BK rank condition not satisfied



  info = 0;
  
  options_ = set_default_option(options_,'loglinear',0);
  options_ = set_default_option(options_,'noprint',0);
  options_ = set_default_option(options_,'olr',0);
  options_ = set_default_option(options_,'olr_beta',1);
  options_ = set_default_option(options_,'qz_criterium',1.000001);

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

tempex = oo_.exo_simul;

it_ = M_.maximum_lag + 1;
if options_.olr
  z = repmat(zeros(M_.endo_nbr,1),1,klen);
else
  z = repmat(dr.ys,1,klen);
end
z = z(iyr0) ;
if options_.order == 1
  [junk,jacobia_] = feval([M_.fname '_dynamic'],z,tempex);
elseif options_.order == 2
    [junk,jacobia_,hessian] = feval([M_.fname '_dynamic'],z,...
				    [oo_.exo_simul ...
		                     oo_.exo_det_simul]);
end

oo_.exo_simul = tempex ;
tempex = [];

% expanding system for Optimal Linear Regulator
if options_.olr
  bet = options_.olr_beta;
  jacobia1 = [];
  n_inst = size(options_.olr_inst,1);

  if ~isfield(olr_state_,'done')
    olr_state_.done = 1;
    olr_state_.old_M_.maximum_lag = M_.maximum_lag;
    olr_state_.old_M_.maximum_lead = M_.maximum_lead;
    olr_state_.old_M_.endo_nbr = M_.endo_nbr;
    olr_state_.old_M_.lead_lag_incidence = M_.lead_lag_incidence;
    
    for i=1:M_.endo_nbr
      temp = ['mult_' int2str(i)];
      lgoo_.endo_simul = strvcat(lgoo_.endo_simul,temp);
    end
    M_.endo_nbr = 2*M_.endo_nbr-n_inst;
    M_.maximum_lag = max(M_.maximum_lag,M_.maximum_lead);
    M_.maximum_lead = M_.maximum_lag;
  end    
  nj = olr_state_.old_M_.endo_nbr-n_inst;
  offset_min = M_.maximum_lag - olr_state_.old_M_.maximum_lag;
  offset_max = M_.maximum_lead - olr_state_.old_M_.maximum_lead;
  newiy = zeros(2*M_.maximum_lag+1,nj+olr_state_.old_M_.endo_nbr);
  jacobia_ = jacobia_(1:nj,:);
  for i=1:2*M_.maximum_lag+1
    if i > offset_min & i <= 2*M_.maximum_lag+1-offset_max
      [junk,k1,k2] = find(olr_state_.old_M_.lead_lag_incidence(i-offset_min,:));
      if i == M_.maximum_lag+1
	jacobia1 = [jacobia1 [jacobia_(:,k2); 2*options_.olr_w]];
      else
	jacobia1 = [jacobia1 [jacobia_(:,k2); ...
		    zeros(olr_state_.old_M_.endo_nbr,length(k1))]];
      end
      newiy(i,k1) = ones(1,length(k1));
    end
    i1  = 2*M_.maximum_lag+2-i;
    if i1 <= 2*M_.maximum_lag+1-offset_max & i1 > offset_min 
      [junk,k1,k2] = find(olr_state_.old_M_.lead_lag_incidence(i1-offset_min,:));
      k3 = find(any(jacobia_(:,k2),2));
      x = zeros(olr_state_.old_M_.endo_nbr,length(k3));
      x(k1,:) = bet^(-i1+M_.maximum_lag+1)*jacobia_(k3,k2)';
      jacobia1  = [jacobia1 [zeros(nj,length(k3)); x]];
      newiy(i,k3+olr_state_.old_M_.endo_nbr) = ones(1,length(k3));
    end      
  end
  jacobia1 = [jacobia1 [jacobia_(:,end-M_.exo_nbr+1:end); ...
		    zeros(olr_state_.old_M_.endo_nbr, M_.exo_nbr)]];
  newiy = newiy';
  newiy = find(newiy(:));
  M_.lead_lag_incidence = zeros(M_.endo_nbr*(M_.maximum_lag+M_.maximum_lead+1),1);
  M_.lead_lag_incidence(newiy) = [1:length(newiy)]';
  M_.lead_lag_incidence =reshape(M_.lead_lag_incidence,M_.endo_nbr,M_.maximum_lag+M_.maximum_lead+1)';
  jacobia_ = jacobia1;
  clear jacobia1
  % computes steady state
  resid = feval([M_.fname '_steady'],zeros(olr_state_.old_M_.endo_nbr,1));
  if resid'*resid < 1e-12
    dr.ys =[dr.ys; zeros(nj,1)];
  else
    AA = zeros(M_.endo_nbr,M_.endo_nbr);
    for i=1:M_.maximum_lag+M_.maximum_lead+1
      [junk,k1,k2] = find(M_.lead_lag_incidence(i,:));
      AA(:,k1) = AA(:,k1)+jacobia_(:,k2);
    end
    dr.ys = -AA\[resid; zeros(nj,1)];
  end
end
% end of code section for Optimal Linear Regulator

klen = M_.maximum_lag + M_.maximum_lead + 1;
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
nz = nnz(M_.lead_lag_incidence);

sdyn = M_.endo_nbr - nstatic;

k0 = M_.lead_lag_incidence(M_.maximum_lag+1,order_var);
k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_lag+1),:);
b = jacobia_(:,k0);

if M_.maximum_lead == 0;  % backward models
  a = jacobia_(:,nonzeros(k1'));
  dr.ghx = zeros(size(a));
  m = 0;
  for i=M_.maximum_lag:-1:1
    k = nonzeros(M_.lead_lag_incidence(i,order_var));
    dr.ghx(:,m+[1:length(k)]) = -b\a(:,k);
    m = m+length(k);
  end
  if M_.exo_nbr
    dr.ghu = -b\jacobia_(:,nz+1:end);
  end
  dr.eigval = eig(transition_matrix(dr));
  dr.rank = 0;
  if any(abs(dr.eigval) > options_.qz_criterium)
    temp = sort(abs(dr.eigval));
    nba = nnz(abs(dr.eigval) > options_.qz_criterium);
    temp = temp(nd-nba+1:nd)-1-options_.qz_criterium;
    info(1) = 3;
    info(2) = temp'*temp;
  end
  return;
end

%forward--looking models
if nstatic > 0
  [Q,R] = qr(b(:,1:nstatic));
  aa = Q'*jacobia_;
else
  aa = jacobia_;
end
a = aa(:,nonzeros(k1'));
b = aa(:,k0);
b10 = b(1:nstatic,1:nstatic);
b11 = b(1:nstatic,nstatic+1:end);
b2 = b(nstatic+1:end,nstatic+1:end);
if any(isinf(a(:)))
  info = 1;
  return
end
if M_.exo_nbr
  fu = aa(:,nz+(1:M_.exo_nbr));
end
clear aa;

% buildind D and E
d = zeros(nd,nd) ;
e = d ;

k = find(kstate(:,2) >= M_.maximum_lag+2 & kstate(:,3));
d(1:sdyn,k) = a(nstatic+1:end,kstate(k,3)) ;
k1 = find(kstate(:,2) == M_.maximum_lag+2);
e(1:sdyn,k1) =  -b2(:,kstate(k1,1)-nstatic);
k = find(kstate(:,2) <= M_.maximum_lag+1 & kstate(:,4));
e(1:sdyn,k) = -a(nstatic+1:end,kstate(k,4)) ;
k2 = find(kstate(:,2) == M_.maximum_lag+1);
k2 = k2(~ismember(kstate(k2,1),kstate(k1,1)));
d(1:sdyn,k2) = b2(:,kstate(k2,1)-nstatic);

if ~isempty(kad)
  for j = 1:size(kad,1)
    d(sdyn+j,kad(j)) = 1 ;
    e(sdyn+j,kae(j)) = 1 ;
  end
end

if ~exist('mjdgges')
  % using Chris Sim's routines
  use_qzdiv = 1;
  [ss,tt,qq,w] = qz(e,d);
  [tt,ss,qq,w] = qzdiv(options_.qz_criterium,tt,ss,qq,w);
  ss1=diag(ss);
  tt1=diag(tt);
  warning_state = warning;
  warning off;
  dr.eigval = ss1./tt1 ;
  warning warning_state;
  nba = nnz(abs(dr.eigval) > options_.qz_criterium);
else
  use_qzdiv = 0;
  [ss,tt,w,sdim,dr.eigval,info1] = mjdgges(e,d,options_.qz_criterium);
  if info1
    info(1) = 2;
    info(2) = info1;
    return
  end
  nba = nd-sdim;
end

nyf = sum(kstate(:,2) > M_.maximum_lag+1);

if task == 1
  dr.rank = rank(w(1:nyf,nd-nyf+1:end));
  dr.eigval = eig(e,d);
  return
end

if nba ~= nyf
  temp = sort(abs(dr.eigval));
  if nba > nyf
    temp = temp(nd-nba+1:nd-nyf)-1-options_.qz_criterium;
    info(1) = 3;
  elseif nba < nyf;
    temp = temp(nd-nyf+1:nd-nba)-1-options_.qz_criterium;
    info(1) = 4;
  end
  info(2) = temp'*temp;
  return
end

np = nd - nyf;
n2 = np + 1;
n3 = nyf;
n4 = n3 + 1;
% derivatives with respect to dynamic state variables
% forward variables
w1 =w(1:n3,n2:nd);
if condest(w1) > 1e9;
  info(1) = 5;
  info(2) = condest(w1);
  return;
else
  gx = -w1'\w(n4:nd,n2:nd)';
end  

% predetermined variables
hx = w(1:n3,1:np)'*gx+w(n4:nd,1:np)';
hx = (tt(1:np,1:np)*hx)\(ss(1:np,1:np)*hx);

k1 = find(kstate(n4:nd,2) == M_.maximum_lag+1);
k2 = find(kstate(1:n3,2) == M_.maximum_lag+2);
dr.ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];

%lead variables actually present in the model
j3 = nonzeros(kstate(:,3));
j4  = find(kstate(:,3));
% derivatives with respect to exogenous variables
if M_.exo_nbr
  a1 = b;
  aa1 = [];
  if nstatic > 0
    aa1 = a1(:,1:nstatic);
  end
  dr.ghu = -[aa1 a(:,j3)*gx(j4,1:npred)+a1(:,nstatic+1:nstatic+ ...
						  npred) a1(:,nstatic+npred+1:end)]\fu;
end

% static variables
if nstatic > 0
  temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
  j5 = find(kstate(n4:nd,4));
  temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
  temp = b10\(temp-b11*dr.ghx);
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

%necessary when using Sims' routines
if use_qzdiv
  gx = real(gx);
  hx = real(hx);
  dr.ghx = real(dr.ghx);
  dr.ghu = real(dr.ghu);
end

%exogenous deterministic variables
if M_.exo_det_nbr > 0
  f1 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_lag+2:end,order_var))));
  f0 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_lag+1,order_var))));
  fudet = sparse(jacobia_(:,nz+M_.exo_nbr+1:end));
  M1 = inv(f0+[zeros(M_.endo_nbr,nstatic) f1*gx zeros(M_.endo_nbr,nyf-nboth)]);
  M2 = M1*f1;
  dr.ghud = cell(M_.exo_det_length,1);
  dr.ghud{1} = -M1*fudet;
  for i = 2:M_.exo_det_length
    dr.ghud{i} = -M2*dr.ghud{i-1}(end-nyf+1:end,:);
  end
end

if options_.order == 1
  return
end

% Second order
%tempex = oo_.exo_simul ;

%hessian = real(hessext('ff1_',[z; oo_.exo_steady_state]))' ;
kk = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+1:end,order_var)),1));
if M_.maximum_lag > 0
  kk = [cumsum(M_.lead_lag_incidence(1:M_.maximum_lag,order_var),1); kk];
end
kk = kk';
kk = find(kk(:));
nk = size(kk,1) + M_.exo_nbr + M_.exo_det_nbr;
k1 = M_.lead_lag_incidence(:,order_var);
k1 = k1';
k1 = k1(:);
k1 = k1(kk);
k2 = find(k1);
kk1(k1(k2)) = k2;
kk1 = [kk1 length(k1)+1:length(k1)+M_.exo_nbr+M_.exo_det_nbr];
kk = reshape([1:nk^2],nk,nk);
kk1 = kk(kk1,kk1);
%[junk,junk,hessian] = feval([M_.fname '_dynamic'],z, oo_.exo_steady_state);
hessian(:,kk1(:)) = hessian;

%oo_.exo_simul = tempex ;
%clear tempex

n1 = 0;
n2 = np;
zx = zeros(np,np);
zu=zeros(np,M_.exo_nbr);
for i=2:M_.maximum_lag+1
  k1 = sum(kstate(:,2) == i);
  zx(n1+1:n1+k1,n2-k1+1:n2)=eye(k1);
  n1 = n1+k1;
  n2 = n2-k1;
end
kk = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+1:end,order_var)),1));
k0 = [1:M_.endo_nbr];
gx1 = dr.ghx;
hu = dr.ghu(nstatic+[1:npred],:);
zx = [zx; gx1];
zu = [zu; dr.ghu];
for i=1:M_.maximum_lead
  k1 = find(kk(i+1,k0) > 0);
  zu = [zu; gx1(k1,1:npred)*hu];
  gx1 = gx1(k1,:)*hx;
  zx = [zx; gx1];
  kk = kk(:,k0);
  k0 = k1;
end
zx=[zx; zeros(M_.exo_nbr,np);zeros(M_.exo_det_nbr,np)];
zu=[zu; eye(M_.exo_nbr);zeros(M_.exo_det_nbr,M_.exo_nbr)];
[n1,n2] = size(zx);
if n1*n1*n2*n2 > 1e7
  rhs = zeros(M_.endo_nbr,n2*n2);
  k1 = 1;
  for i1 = 1:n2
      for i2 = 1:n2
	rhs(:,k1) = hessian*kron(zx(:,i1),zx(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zx,zx);
end
rhs = -rhs;

%lhs
n = M_.endo_nbr+sum(kstate(:,2) > M_.maximum_lag+1 & kstate(:,2) < M_.maximum_lag+M_.maximum_lead+1);
A = zeros(n,n);
B = zeros(n,n);
A(1:M_.endo_nbr,1:M_.endo_nbr) = jacobia_(:,M_.lead_lag_incidence(M_.maximum_lag+1,order_var));
% variables with the highest lead
k1 = find(kstate(:,2) == M_.maximum_lag+M_.maximum_lead+1);
if M_.maximum_lead > 1
  k2 = find(kstate(:,2) == M_.maximum_lag+M_.maximum_lead);
  [junk,junk,k3] = intersect(kstate(k1,1),kstate(k2,1));
else
  k2 = [1:M_.endo_nbr];
  k3 = kstate(k1,1);
end
% Jacobian with respect to the variables with the highest lead
B(1:M_.endo_nbr,end-length(k2)+k3) = jacobia_(:,kstate(k1,3)+M_.endo_nbr);
offset = M_.endo_nbr;
k0 = [1:M_.endo_nbr];
gx1 = dr.ghx;
for i=1:M_.maximum_lead-1
  k1 = find(kstate(:,2) == M_.maximum_lag+i+1);
  [k2,junk,k3] = find(kstate(k1,3));
  A(1:M_.endo_nbr,offset+k2) = jacobia_(:,k3+M_.endo_nbr);
  n1 = length(k1);
  A(offset+[1:n1],nstatic+[1:npred]) = -gx1(kstate(k1,1),1:npred);
  gx1 = gx1*hx;
  A(offset+[1:n1],offset+[1:n1]) = eye(n1);
  n0 = length(k0);
  E = eye(n0);
  if i == 1
    [junk,junk,k4]=intersect(kstate(k1,1),[1:M_.endo_nbr]);
  else
    [junk,junk,k4]=intersect(kstate(k1,1),kstate(k0,1));
  end
  i1 = offset-n0+n1;
  B(offset+[1:n1],offset-n0+[1:n0]) = -E(k4,:);
  k0 = k1;
  offset = offset + n1;
end
[junk,k1,k2] = find(M_.lead_lag_incidence(M_.maximum_lag+M_.maximum_lead+1,order_var));
A(1:M_.endo_nbr,nstatic+1:nstatic+npred)=...
    A(1:M_.endo_nbr,nstatic+[1:npred])+jacobia_(:,k2)*gx1(k1,1:npred);
C = hx;
D = [rhs; zeros(n-M_.endo_nbr,size(rhs,2))];
dr.ghxx = gensylv(2,A,B,C,D);
if exist('gensylv')
  dr.ghxx = gensylv(2,A,B,C,D);
else
  C = kron(hx,hx); 
  x0 = sylvester3(A,B,C,D);
  dr.ghxx = sylvester3a(x0,A,B,C,D);
end

%ghxu
%rhs
hu = dr.ghu(nstatic+1:nstatic+npred,:);
%kk = reshape([1:np*np],np,np);
%kk = kk(1:npred,1:npred);
%rhs = -hessian*kron(zx,zu)-f1*dr.ghxx(end-nyf+1:end,kk(:))*kron(hx(1:npred,:),hu(1:npred,:));
if n1*n1*n2*M_.exo_nbr > 1e7
  rhs = zeros(M_.endo_nbr,n2*M_.exo_nbr);
  k1 = 1;
  for i1 = 1:n2
      for i2 = 1:M_.exo_nbr
	rhs(:,k1) = hessian*kron(zx(:,i1),zu(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zx,zu);
end
nyf1 = sum(kstate(:,2) == M_.maximum_lag+2);
hu1 = [hu;zeros(np-npred,M_.exo_nbr)];
%B1 = [B(1:M_.endo_nbr,:);zeros(size(A,1)-M_.endo_nbr,size(B,2))];
rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B*dr.ghxx*kron(hx,hu1);


%lhs
dr.ghxu = A\rhs;

%ghuu
%rhs
kk = reshape([1:np*np],np,np);
kk = kk(1:npred,1:npred);
if n1*n1*M_.exo_nbr*M_.exo_nbr > 1e7
  rhs = zeros(M_.endo_nbr,M_.exo_nbr*M_.exo_nbr);
  k1 = 1;
  for i1 = 1:M_.exo_nbr
      for i2 = 1:M_.exo_nbr
	rhs(:,k1) = hessian*kron(zu(:,i1),zu(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zu,zu);
end

rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B*dr.ghxx*kron(hu1,hu1);

%lhs
dr.ghuu = A\rhs;

dr.ghxx = dr.ghxx(1:M_.endo_nbr,:);
dr.ghxu = dr.ghxu(1:M_.endo_nbr,:);
dr.ghuu = dr.ghuu(1:M_.endo_nbr,:);


% dr.ghs2
% derivatives of F with respect to forward variables
% reordering predetermined variables in diminishing lag order
O1 = zeros(M_.endo_nbr,nstatic);
O2 = zeros(M_.endo_nbr,M_.endo_nbr-nstatic-npred);
LHS = jacobia_(:,M_.lead_lag_incidence(M_.maximum_lag+1,order_var));
RHS = zeros(M_.endo_nbr,M_.exo_nbr^2);
kk = find(kstate(:,2) == M_.maximum_lag+2);
gu = dr.ghu; 
guu = dr.ghuu; 
Gu = [dr.ghu(nstatic+[1:npred],:); zeros(np-npred,M_.exo_nbr)];
Guu = [dr.ghuu(nstatic+[1:npred],:); zeros(np-npred,M_.exo_nbr*M_.exo_nbr)];
E = eye(M_.endo_nbr);
M_.lead_lag_incidenceordered = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_lag+1:end,order_var)),1));
if M_.maximum_lag > 0
  M_.lead_lag_incidenceordered = [cumsum(M_.lead_lag_incidence(1:M_.maximum_lag,order_var),1); M_.lead_lag_incidenceordered];
end
M_.lead_lag_incidenceordered = M_.lead_lag_incidenceordered';
M_.lead_lag_incidenceordered = M_.lead_lag_incidenceordered(:);
k1 = find(M_.lead_lag_incidenceordered);
M_.lead_lag_incidenceordered(k1) = [1:length(k1)]';
M_.lead_lag_incidenceordered =reshape(M_.lead_lag_incidenceordered,M_.endo_nbr,M_.maximum_lag+M_.maximum_lead+1)';
kh = reshape([1:nk^2],nk,nk);
kp = sum(kstate(:,2) <= M_.maximum_lag+1);
E1 = [eye(npred); zeros(kp-npred,npred)];
H = E1;
hxx = dr.ghxx(nstatic+[1:npred],:);
for i=1:M_.maximum_lead
  for j=i:M_.maximum_lead
    [junk,k2a,k2] = find(M_.lead_lag_incidence(M_.maximum_lag+j+1,order_var));
    [junk,k3a,k3] = find(M_.lead_lag_incidenceordered(M_.maximum_lag+j+1,:));
    RHS = RHS + jacobia_(:,k2)*guu(k2a,:)+hessian(:,kh(k3,k3))* ...
	  kron(gu(k3a,:),gu(k3a,:));
  end

  % LHS
  [junk,k2a,k2] = find(M_.lead_lag_incidence(M_.maximum_lag+i+1,order_var));
  LHS = LHS + jacobia_(:,k2)*(E(k2a,:)+[O1(k2a,:) dr.ghx(k2a,:)*H O2(k2a,:)]);
  
  if i == M_.maximum_lead 
    break
  end
  
  kk = find(kstate(:,2) == M_.maximum_lag+i+1);
  gu = dr.ghx*Gu;
  GuGu = kron(Gu,Gu);
  guu = dr.ghx*Guu+dr.ghxx*GuGu;
  Gu = hx*Gu;
  Guu = hx*Guu;
  Guu(end-npred+1:end,:) = Guu(end-npred+1:end,:) + hxx*GuGu;

  H = E1 + hx*H;
end
RHS = RHS*M_.Sigma_e(:);
dr.fuu = RHS;
RHS = -RHS-dr.fbias;
dr.ghs2 = LHS\RHS;

% deterministic exogenous variables
if M_.exo_det_nbr > 0
  hud = dr.ghud{1}(nstatic+1:nstatic+npred,:);
  zud=[zeros(np,M_.exo_det_nbr);dr.ghud{1};gx(:,1:npred)*hud;zeros(M_.exo_nbr,M_.exo_det_nbr);eye(M_.exo_det_nbr)];
  R1 = hessian*kron(zx,zud);
  dr.ghxud = cell(M_.exo_det_length,1);
  kf = [M_.endo_nbr-nyf+1:M_.endo_nbr];
  kp = nstatic+[1:npred];
  dr.ghxud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{1}(kp,:)));
  Eud = eye(M_.exo_det_nbr);
  for i = 2:M_.exo_det_length
    hudi = dr.ghud{i}(kp,:);
    zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
    R2 = hessian*kron(zx,zudi);
    dr.ghxud{i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hx,Eud)+dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{i}(kp,:)))-M1*R2;
  end
  R1 = hessian*kron(zu,zud);
  dr.ghudud = cell(M_.exo_det_length,1);
  kf = [M_.endo_nbr-nyf+1:M_.endo_nbr];

  dr.ghuud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghu(kp,:),dr.ghud{1}(kp,:)));
  Eud = eye(M_.exo_det_nbr);
  for i = 2:M_.exo_det_length
    hudi = dr.ghud{i}(kp,:);
    zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
    R2 = hessian*kron(zu,zudi);
    dr.ghuud{i} = -M2*dr.ghxud{i-1}(kf,:)*kron(hu,Eud)-M1*R2;
  end
  R1 = hessian*kron(zud,zud);
  dr.ghudud = cell(M_.exo_det_length,M_.exo_det_length);
  dr.ghudud{1,1} = -M1*R1-M2*dr.ghxx(kf,:)*kron(hud,hud);
  for i = 2:M_.exo_det_length
    hudi = dr.ghud{i}(nstatic+1:nstatic+npred,:);
    zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi+dr.ghud{i-1}(kf,:);zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
    R2 = hessian*kron(zudi,zudi);
    dr.ghudud{i,i} = -M2*(dr.ghudud{i-1,i-1}(kf,:)+...
			  2*dr.ghxud{i-1}(kf,:)*kron(hudi,Eud) ...
			  +dr.ghxx(kf,:)*kron(hudi,hudi))-M1*R2;
    R2 = hessian*kron(zud,zudi);
    dr.ghudud{1,i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hud,Eud)+...
			  dr.ghxx(kf,:)*kron(hud,hudi))...
	-M1*R2;
    for j=2:i-1
      hudj = dr.ghud{j}(kp,:);
      zudj=[zeros(np,M_.exo_det_nbr);dr.ghud{j};gx(:,1:npred)*hudj;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
      R2 = hessian*kron(zudj,zudi);
      dr.ghudud{j,i} = -M2*(dr.ghudud{j-1,i-1}(kf,:)+dr.ghxud{j-1}(kf,:)* ...
			    kron(hudi,Eud)+dr.ghxud{i-1}(kf,:)* ...
			    kron(hudj,Eud)+dr.ghxx(kf,:)*kron(hudj,hudi))-M1*R2;
    end
    
  end
end
% 01/08/2001 MJ put return on iorder == 1 after defining dr.kstate and dr.kdyn
% 01/17/2001 MJ added dr.delta_s: correction factor for order = 2
% 01/21/2001 FC correction of delta_s for more than 1 shock
% 01/23/2001 MJ suppression of redundant sum() in delta_s formula
% 02/22/2001 MJ stderr_ replaced by Sigma_e_
% 08/24/2001 MJ changed the order of the variables, separates static
%               variables and handles only one instance of variables both
%               in lead and lag
% 08/24/2001 MJ added sigma to state variables as in Schmitt-Grohe and
%               Uribe (2001)
% 10/20/2002 MJ corrected lags on several periods bug
% 10/30/2002 MJ corrected lags on several periods bug on static when some
%               intermediary lags are missing
% 12/08/2002 MJ uses sylvester3 to solve for dr.ghxx
% 01/01/2003 MJ added dr.fbias for iterative for dr_algo == 1
% 02/21/2003 MJ corrected bug for models without lagged variables
% 03/02/2003 MJ fixed second order for lag on several periods
% 05/21/2003 MJ add check call argument and make computation for CHECK
% 06/01/2003 MJ added a test for M_.maximum_lead > 1 and order > 1
% 08/28/2003 MJ corrected bug in computation of 2nd order (ordering of
%               forward variable in LHS)
% 08/29/2003 MJ use Sims routine if mjdgges isn't available

   





