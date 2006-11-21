function [Vi, Vk, Vkk] =vardec(k, varargin);
global oo_  M_   

lgx_orig_ord_=M_.exo_names_orig_ord;
Sigma_e_=M_.Sigma_e;
lgy_=M_.endo_names;
lgx_=M_.exo_names;
fname_=M_.fname;
dr_=oo_.dr;
if isempty(k), k=40; end,
nvar = length(varargin);
if nvar == 0
    nvar  = length(dr_.order_var);
    ivar  = transpose(1:nvar);
    var_list = lgy_;
else
    ivar=zeros(nvar,1);
    for i=1:nvar
        if i>1,
            var_list=str2mat(var_list,varargin{i});
        else
            var_list(i,:)=varargin{i};
        end
        i_tmp = strmatch(var_list(i,:),lgy_,'exact');
        if isempty(i_tmp)
            error (['One of the variable specified does not exist']) ;
        else
            ivar(i) = i_tmp;
        end	
    end
end

[T,B,SteadyState,info] = dynare_resolve;

[dum, idum]=sort(dr_.order_var);
%BB=B(idum,:);
S0=diag(diag(Sigma_e_));
V=0.*T;
Vi=zeros(size(lgy_,1), size(lgx_,1), k);
Vii=zeros(size(T,1), size(T,2),size(lgx_,1));
%T=dr_.ghx;
%B=dr_.ghu;
for i=1:k,
    V=T*V*T'+B*Sigma_e_*B';
    V0=diag(V);
    V0=V0(idum);
    Vdum = V(idum, idum);
    Vk(:,i)=V0(ivar);
    Vkk(:,:,i)=Vdum(ivar, ivar);
    iv{i}=find(V0(ivar)>1.e-12);
    iv0=find(V0<=1.e-12);
    V0(iv0)=ones(length(iv0),1);
    for j=1:length(Sigma_e_),    
        Vii(:,:,j)= T*Vii(:,:,j)*T' + B(:,j)*S0(j,j)*B(:,j)';
        V0i=diag(Vii(:,:,j));
        Vi(:,j,i)=V0i(idum)./V0;
    end
end

BQ=B*Sigma_e_*B';
Vinf = lyapunov_symm(T,BQ);
V0=diag(Vinf);
V0=V0(idum);
ivinf=find(V0(ivar)>1.e-12);
iv0=find(V0<=1.e-12);
V0(iv0)=ones(length(iv0),1);

for j=1:length(Sigma_e_),   
    BQ=B(:,j)*S0(j,j)*B(:,j)';
    V00 = lyapunov_symm(T,BQ);
    V00  = diag(V00);
    Vinfi(:,j) = V00(idum)./V0;
    %Vinfi(:,j) = inv(eye(length(T))-kron(T,T))*BQ(:)./V0;
end

iprint=[1, 4, 8, 16, 40, 100, 200, 400];
lh = size(deblank(lgy_(ivar,:)),2)+2;
headers = lgx_;
%headers(lgx_orig_ord_,:) = headers;
headers = strvcat(' ',headers);
for j=1:length(iprint),
    
    ip=iprint(j);
    if k>=ip,
        title=['VARIANCE DECOMPOSITION (',int2str(ip),'-step ahead, in percent)'];
        
        
        table(title,headers,deblank(lgy_(ivar(iv{ip}),:)),100*Vi(ivar(iv{ip}),:,ip), ...
            lh,8,2);
    end
end

title=['VARIANCE DECOMPOSITION (inf-step ahead, in percent)'];
        
table(title,headers,deblank(lgy_(ivar(ivinf),:)),100*Vinfi(ivar(ivinf),:), ...
    lh,8,2);

save([fname_,'_vardec'],'Vi', 'Vinfi') 
