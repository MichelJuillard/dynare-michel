function [H, dA, dOm, Hss, gp, d2A, d2Om, H2ss] = getH(A, B, M_,oo_,options_,kronflag,indx,indexo)

% computes derivative of reduced form linear model w.r.t. deep params
%
% Copyright (C) 2010-2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if nargin<3 || isempty(kronflag), kronflag = 0; end
if nargin<4 || isempty(indx), indx = [1:M_.param_nbr]; end,
if nargin<5 || isempty(indexo), indexo = []; end,

[I,J]=find(M_.lead_lag_incidence');
yy0=oo_.dr.ys(I);
param_nbr = length(indx);
if nargout>5,
    param_nbr_2 = param_nbr*(param_nbr+1)/2;
end

m = size(A,1);
n = size(B,2);

if kronflag==-2,
    if nargout>5,
        [residual, g1, g2 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
            M_.params, oo_.dr.ys, 1);
        g22 = hessian('thet2tau',[M_.params(indx)],options_.gstep,M_, oo_, indx,[],-1);
        H2ss=g22(1:M_.endo_nbr,:);
        H2ss = reshape(H2ss,[M_.endo_nbr param_nbr param_nbr]);
        g22=g22(M_.endo_nbr+1:end,:);
        g22 = reshape(g22,[size(g1) param_nbr param_nbr]);
    else
        [residual, g1 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
            M_.params, oo_.dr.ys, 1);        
    end
    gp = fjaco('thet2tau',[M_.params(indx)],M_, oo_, indx,[],-1);
    Hss=gp(1:M_.endo_nbr,:);
    gp=gp(M_.endo_nbr+1:end,:);
    gp = reshape(gp,[size(g1) param_nbr]);
else

% yy0=[];
% for j=1:size(M_.lead_lag_incidence,1);
%     yy0 = [ yy0; oo_.dr.ys(find(M_.lead_lag_incidence(j,:)))];
% end
dyssdtheta=zeros(length(oo_.dr.ys),M_.param_nbr);
d2yssdtheta=zeros(length(oo_.dr.ys),M_.param_nbr,M_.param_nbr);
if nargout>5,
    [residual, gg1, gg2] = feval([M_.fname,'_static'],oo_.dr.ys, oo_.exo_steady_state', M_.params);
    [df, gp, d2f] = feval([M_.fname,'_params_derivs'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1, dyssdtheta, d2yssdtheta);
    d2f = get_all_resid_2nd_derivs(d2f,length(oo_.dr.ys),M_.param_nbr);
    d2f = d2f(:,indx,indx);
    if isempty(find(gg2)),
        for j=1:length(indx),
        d2yssdtheta(:,indx,j) = -gg1\d2f(:,:,j);
        end
    else
        gam = d2f*0;
        for j=1:length(indx),
        d2yssdtheta(:,indx,j) = -gg1\(d2f(:,:,j)+gam(:,:,j));
        end
    end
else
    [residual, gg1] = feval([M_.fname,'_static'],oo_.dr.ys, oo_.exo_steady_state', M_.params);
    df = feval([M_.fname,'_params_derivs'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1, dyssdtheta, d2yssdtheta);
end
dyssdtheta = -gg1\df;

if any(any(isnan(dyssdtheta))),    
    [U,T] = schur(gg1);
    qz_criterium=options_.qz_criterium;
    e1 = abs(ordeig(T)) < qz_criterium-1;
    k = sum(e1);       % Number non stationary variables.
%     n = length(e1)-k;  % Number of stationary variables.
    [U,T] = ordschur(U,T,e1);
    T = T(k+1:end,k+1:end);
    dyssdtheta = -U(:,k+1:end)*(T\U(:,k+1:end)')*df;
    if nargout>5,
        for j=1:length(indx),
            d2yssdtheta(:,:,j) = -U(:,k+1:end)*(T\U(:,k+1:end)')*d2f(:,:,j);
        end
    end
end
if nargout>5,
    [df, gp, d2f, gpp] = feval([M_.fname,'_params_derivs'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1, dyssdtheta, d2yssdtheta);
    H2ss = d2yssdtheta(oo_.dr.order_var,:,:);
    [residual, g1, g2, g3] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1);
    nelem=size(g1,2);
    g22 = get_all_2nd_derivs(gpp,m,nelem,M_.param_nbr);
    g22 = g22(:,:,indx,indx);
else
    [df, gp] = feval([M_.fname,'_params_derivs'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1, dyssdtheta,d2yssdtheta);
    [residual, g1, g2 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', ...
        M_.params, oo_.dr.ys, 1);
end

Hss = dyssdtheta(oo_.dr.order_var,indx);
dyssdtheta = dyssdtheta(I,:);
[nr, nc]=size(g2);
nc = sqrt(nc);
ns = max(max(M_.lead_lag_incidence));
gp2 = gp*0;
for j=1:nr,
    [I J]=ind2sub([nc nc],find(g2(j,:)));
    for i=1:nc,
        is = find(I==i);
        is = is(find(J(is)<=ns));
        if ~isempty(is),
            g20=full(g2(j,find(g2(j,:))));
            gp2(j,i,:)=g20(is)*dyssdtheta(J(is),:);
        end
    end
end


gp = gp+gp2;
gp = gp(:,:,indx);
end



klen = M_.maximum_endo_lag + M_.maximum_endo_lead + 1;
k11 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_endo_lag+1),:);
a = g1(:,nonzeros(k11'));
da = gp(:,nonzeros(k11'),:);
if nargout > 5,
    d2a = g22(:,nonzeros(k11'),:,:);
end
kstate = oo_.dr.kstate;

GAM1 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg1 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);
% nf = find(M_.lead_lag_incidence(M_.maximum_endo_lag+2,:));
% GAM1(:, nf) = -g1(:,M_.lead_lag_incidence(M_.maximum_endo_lag+2,nf));

k1 = find(kstate(:,2) == M_.maximum_endo_lag+2 & kstate(:,3));
GAM1(:, kstate(k1,1)) = -a(:,kstate(k1,3));
Dg1(:, kstate(k1,1), :) = -da(:,kstate(k1,3),:);
if nargout > 5,
    D2g1 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr,param_nbr);
    D2g1(:, kstate(k1,1), :, :) = -d2a(:,kstate(k1,3),:,:);
end
% k = find(kstate(:,2) > M_.maximum_endo_lag+2 & kstate(:,3));
% kk = find(kstate(:,2) > M_.maximum_endo_lag+2 & ~kstate(:,3));
% nauxe = 0;
% if ~isempty(k),
%     ax(:,k)= -a(:,kstate(k,3));
%     ax(:,kk)= 0;
%     dax(:,k,:)= -da(:,kstate(k,3),:);
%     dax(:,kk,:)= 0;
%     nauxe=size(ax,2);
%     GAM1 = [ax GAM1];
%     Dg1 = cat(2,dax,Dg1);
%     clear ax
% end


[junk,cols_b,cols_j] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, ...
                                                  oo_.dr.order_var));
GAM0 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg0 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);
GAM0(:,cols_b) = g1(:,cols_j);
Dg0(:,cols_b,:) = gp(:,cols_j,:);
if nargout > 5,
    D2g0 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr,param_nbr);
    D2g0(:, cols_b, :, :) = g22(:,cols_j,:,:);
end
% GAM0 = g1(:,M_.lead_lag_incidence(M_.maximum_endo_lag+1,:));


k2 = find(kstate(:,2) == M_.maximum_endo_lag+1 & kstate(:,4));
GAM2 = zeros(M_.endo_nbr,M_.endo_nbr);
Dg2 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr);
% nb = find(M_.lead_lag_incidence(1,:));
% GAM2(:, nb) = -g1(:,M_.lead_lag_incidence(1,nb));
GAM2(:, kstate(k2,1)) = -a(:,kstate(k2,4));
Dg2(:, kstate(k2,1), :) = -da(:,kstate(k2,4),:);
if nargout > 5,
    D2g2 = zeros(M_.endo_nbr,M_.endo_nbr,param_nbr,param_nbr);
    D2g2(:, kstate(k2,1), :, :) = -d2a(:,kstate(k2,4),:,:);
end% naux = 0;
% k = find(kstate(:,2) < M_.maximum_endo_lag+1 & kstate(:,4));
% kk = find(kstate(:,2) < M_.maximum_endo_lag+1 );
% if ~isempty(k),
%     ax(:,k)= -a(:,kstate(k,4));
%     ax = ax(:,kk);
%     dax(:,k,:)= -da(:,kstate(k,4),:);
%     dax = dax(:,kk,:);
%     naux = size(ax,2);
%     GAM2 = [GAM2 ax];
%     Dg2 = cat(2,Dg2,dax);
% end
% 
% GAM0 = blkdiag(GAM0,eye(naux));
% Dg0 = cat(1,Dg0,zeros(naux+nauxe,M_.endo_nbr,param_nbr));
% Dg0 = cat(2,Dg0,zeros(m+nauxe,naux,param_nbr));
% Dg0 = cat(2,zeros(m+nauxe,nauxe,param_nbr),Dg0);
% 
% GAM2 = [GAM2 ; A(M_.endo_nbr+1:end,:)];
% Dg2 = cat(1,Dg2,zeros(naux+nauxe,m,param_nbr));
% Dg2 = cat(2,zeros(m+nauxe,nauxe,param_nbr),Dg2);
% GAM2 = [zeros(m+nauxe,nauxe) [GAM2; zeros(nauxe,m)]];
% GAM0 = [[zeros(m,nauxe);(eye(nauxe))] [GAM0; zeros(nauxe,m)]];

GAM3 = -g1(:,length(yy0)+1:end);
% GAM3 = -g1(oo_.dr.order_var,length(yy0)+1:end);
% GAM3 = [GAM3; zeros(naux+nauxe,size(GAM3,2))];
% Dg3 = -gp(oo_.dr.order_var,length(yy0)+1:end,:);
Dg3 = -gp(:,length(yy0)+1:end,:);
if nargout>5,
    D2g3 = -g22(:,length(yy0)+1:end,:,:);
    clear g22 d2a
end
% Dg3 = cat(1,Dg3,zeros(naux+nauxe,size(GAM3,2),param_nbr));

% auxe = zeros(0,1);
% k0 = kstate(find(kstate(:,2) >= M_.maximum_endo_lag+2),:);;
% i0 = find(k0(:,2) == M_.maximum_endo_lag+2);
% for i=M_.maximum_endo_lag+3:M_.maximum_endo_lag+2+M_.maximum_endo_lead,
%     i1 = find(k0(:,2) == i);
%     n1 = size(i1,1);
%     j = zeros(n1,1);
%     for j1 = 1:n1
%         j(j1) = find(k0(i0,1)==k0(i1(j1),1));
%     end
%     auxe = [auxe; i0(j)];
%     i0 = i1;
% end
% auxe = [(1:size(auxe,1))' auxe(end:-1:1)];
% n_ir1 = size(auxe,1);
% nr = m + n_ir1;
% GAM1 = [[GAM1 zeros(size(GAM1,1),naux)]; zeros(naux+nauxe,m+nauxe)];
% GAM1(m+1:end,:) = sparse(auxe(:,1),auxe(:,2),ones(n_ir1,1),n_ir1,nr);
% Dg1 = cat(2,Dg1,zeros(M_.endo_nbr,naux,param_nbr));
% Dg1 = cat(1,Dg1,zeros(naux+nauxe,m+nauxe,param_nbr));

% Ax = A;
% A1 = A;
% Bx = B;
% B1 = B;
% for j=1:M_.maximum_endo_lead-1,
%     A1 = A1*A;
%     B1 = A*B1;
%     k = find(kstate(:,2) == M_.maximum_endo_lag+2+j );
%     Ax = [A1(kstate(k,1),:); Ax];
%     Bx = [B1(kstate(k,1),:); Bx];
% end
% Ax = [zeros(m+nauxe,nauxe) Ax];
% A0 = A;
% A=Ax; clear Ax A1;
% B0=B;
% B = Bx; clear Bx B1;

% m = size(A,1);
% n = size(B,2);

% Dg1 = zeros(m,m,param_nbr);
% Dg1(:, nf, :) = -gp(:,M_.lead_lag_incidence(3,nf),:);

% Dg0 = gp(:,M_.lead_lag_incidence(2,:),:);

% Dg2 = zeros(m,m,param_nbr);
% Dg2(:, nb, :) = -gp(:,M_.lead_lag_incidence(1,nb),:);

% Dg3 = -gp(:,length(yy0)+1:end,:);

if kronflag==1, % kronecker products
    Dg0=reshape(Dg0,m^2,param_nbr);
    Dg1=reshape(Dg1,m^2,param_nbr);
    Dg2=reshape(Dg2,m^2,param_nbr);
    for j=1:param_nbr,
        Dg3(:,:,j)=Dg3(:,:,j)*M_.Sigma_e;
    end
    Dg3=reshape(Dg3,m*n,param_nbr);
    Om = B*M_.Sigma_e*B';
    Im = eye(m);
    Dm = duplication(m);
    DmPl = inv(Dm'*Dm)*Dm';
    Kmm = commutation(m,m);
    Kmn = commutation(m,n);


    Da = [eye(m^2),zeros(m^2,m*(m+1)/2)];
    Dom = [zeros(m*(m+1)/2,m^2),eye(m*(m+1)/2)];


    Df1Dtau = ( kron(Im,GAM0) - kron(A',GAM1) - kron(Im,GAM1*A) )*Da;

    Df1Dthet = kron(A',Im)*Dg0 - kron( (A')^2,Im)*Dg1 - Dg2;

    Df2Dtau = DmPl*( kron(GAM0,GAM0) - kron(GAM0,GAM1*A) - kron(GAM1*A,GAM0) + kron(GAM1*A,GAM1*A) )*Dm*Dom - ...
              DmPl*( kron(GAM0*Om,GAM1) + kron(GAM1,GAM0*Om)*Kmm - kron(GAM1*A*Om,GAM1) - kron(GAM1,GAM1*A*Om)*Kmm )*Da;


    Df2Dthet = DmPl*( kron(GAM0*Om,Im) + kron(Im,GAM0*Om)*Kmm - kron(Im,GAM1*A*Om)*Kmm - kron(GAM1*A*Om,Im) )*Dg0 - ...
        DmPl*( kron(GAM0*Om*A',Im) + kron(Im,GAM0*Om*A')*Kmm - kron(Im,GAM1*A*Om*A')*Kmm - kron(GAM1*A*Om*A',Im) )*Dg1 -...
        DmPl*( kron(GAM3,Im) + kron(Im,GAM3)*Kmn )*Dg3;


    DfDtau  = [Df1Dtau;Df2Dtau];
    DfDthet = [Df1Dthet;Df2Dthet];

    H = -DfDtau\DfDthet;
    x = reshape(H(1:m*m,:),m,m,param_nbr);
    y = reshape(Dm*H(m*m+1:end,:),m,m,param_nbr);
    if nauxe,
        x = x(nauxe+1:end,nauxe+1:end,:);
        y = y(nauxe+1:end,nauxe+1:end,:);
        dA = x;
        dOm = y;
        m = size(y,1);
        x = reshape(x,m*m,param_nbr);
        Dm = duplication(m);
        DmPl = inv(Dm'*Dm)*Dm';
        y = DmPl*reshape(y,m*m,param_nbr);
        H = [x;y];
    else
        dA = x;
        dOm = y;
    end

    Hx = [];
    if ~isempty(indexo),
        dSig = zeros(M_.exo_nbr,M_.exo_nbr);
        dOm = cat(3,zeros(size(dOm,1),size(dOm,1),length(indexo)),dOm);
        for j=1:length(indexo)
            dSig(indexo(j),indexo(j)) = 2*sqrt(M_.Sigma_e(indexo(j),indexo(j)));
            y = B*dSig*B';
            y = y(nauxe+1:end,nauxe+1:end);
            Hx(:,j) = [zeros((m-nauxe)^2,1); dyn_vech(y)];
            if nargout>1,
                dOm(:,:,j) = y;
            end
            dSig(indexo(j),indexo(j)) = 0;
        end
    end
    H = [ [zeros(M_.endo_nbr,length(indexo)) Hss]; [Hx H]];

elseif kronflag==-1, % perturbation
    fun = 'thet2tau';
    params0 = M_.params;
    H = fjaco(fun,[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)], M_, oo_, indx, indexo);
    assignin('base','M_', M_);
    assignin('base','oo_', oo_);

else % generalized sylvester equation

    % solves a*x+b*x*c=d
    a = (GAM0-GAM1*A);
    inva = inv(a);
    b = -GAM1;
    c = A;
    elem = zeros(m,m,param_nbr);
    d = elem;
    for j=1:param_nbr,
        elem(:,:,j) = (Dg0(:,:,j)-Dg1(:,:,j)*A);
        d(:,:,j) = Dg2(:,:,j)-elem(:,:,j)*A;
    end
    xx=sylvester3(a,b,c,d);
    flag=1;
    icount=0;
    while flag && icount<4,
        [xx, flag]=sylvester3a(xx,a,b,c,d);
        icount=icount+1;
    end
    H=zeros(m*m+m*(m+1)/2,param_nbr+length(indexo));
    if nargout>1,
        dOm = zeros(m,m,param_nbr+length(indexo));
        dA=zeros(m,m,param_nbr+length(indexo));
        dB=zeros(m,n,param_nbr);
    end
    if ~isempty(indexo),
        dSig = zeros(M_.exo_nbr,M_.exo_nbr,length(indexo));
        for j=1:length(indexo)
            dSig(indexo(j),indexo(j),j) = 2*sqrt(M_.Sigma_e(indexo(j),indexo(j)));
            y = B*dSig(:,:,j)*B';
%             y = y(nauxe+1:end,nauxe+1:end);
%             H(:,j) = [zeros((m-nauxe)^2,1); dyn_vech(y)];
            H(:,j) = [zeros(m^2,1); dyn_vech(y)];
            if nargout>1,
                dOm(:,:,j) = y;
            end
%             dSig(indexo(j),indexo(j)) = 0;
        end
    end
    for j=1:param_nbr,
        x = xx(:,:,j);
        y = inva * (Dg3(:,:,j)-(elem(:,:,j)-GAM1*x)*B);
        if nargout>1,
            dB(:,:,j) = y;
        end
        y = y*M_.Sigma_e*B'+B*M_.Sigma_e*y';
%         x = x(nauxe+1:end,nauxe+1:end);
%         y = y(nauxe+1:end,nauxe+1:end);
        if nargout>1,
            dA(:,:,j+length(indexo)) = x;
            dOm(:,:,j+length(indexo)) = y;
        end
        H(:,j+length(indexo)) = [x(:); dyn_vech(y)];
    end
    %     for j=1:param_nbr,
    %         disp(['Derivatives w.r.t. ',M_.param_names(indx(j),:),', ',int2str(j),'/',int2str(param_nbr)])
    %         elem = (Dg0(:,:,j)-Dg1(:,:,j)*A);
    %         d = Dg2(:,:,j)-elem*A;
    %         x=sylvester3(a,b,c,d);
    % %         x=sylvester3a(x,a,b,c,d);
    %         y = inva * (Dg3(:,:,j)-(elem-GAM1*x)*B);
    %         y = y*B'+B*y';
    %         x = x(nauxe+1:end,nauxe+1:end);
    %         y = y(nauxe+1:end,nauxe+1:end);
    %         H(:,j) = [x(:); dyn_vech(y)];
    %     end
    H = [[zeros(M_.endo_nbr,length(indexo)) Hss]; H];

end

if nargout > 5,
    tot_param_nbr = param_nbr + length(indexo);
    tot_param_nbr_2 = tot_param_nbr*(tot_param_nbr+1)/2;
    d = zeros(m,m,param_nbr_2);
    d2A = zeros(m,m,tot_param_nbr,tot_param_nbr);
    d2Om = zeros(m,m,tot_param_nbr,tot_param_nbr);
    d2B = zeros(m,n,tot_param_nbr,tot_param_nbr);
    cc=triu(ones(param_nbr,param_nbr));
    [i2,j2]=find(cc);
    cc = blkdiag( zeros(length(indexo),length(indexo)), cc);
    [ipar2]=find(cc);
    ctot=triu(ones(tot_param_nbr,tot_param_nbr));
    ctot(1:length(indexo),1:length(indexo))=eye(length(indexo));
    [itot2, jtot2]=find(ctot);
    jcount=0;
    for j=1:param_nbr,
        for i=j:param_nbr,
        elem1 = (D2g0(:,:,j,i)-D2g1(:,:,j,i)*A);
        elem1 = D2g2(:,:,j,i)-elem1*A;
        elemj0 = Dg0(:,:,j)-Dg1(:,:,j)*A;
        elemi0 = Dg0(:,:,i)-Dg1(:,:,i)*A;
        elem2 = -elemj0*xx(:,:,i)-elemi0*xx(:,:,j);
        elem2 = elem2 + ( Dg1(:,:,j)*xx(:,:,i) + Dg1(:,:,i)*xx(:,:,j) )*A;
        elem2 = elem2 + GAM1*( xx(:,:,i)*xx(:,:,j) + xx(:,:,j)*xx(:,:,i));
        jcount=jcount+1;
        d(:,:,jcount) = elem1+elem2;
        end
    end
    xx2=sylvester3(a,b,c,d);
    flag=1;
    icount=0;
    while flag && icount<4,
        [xx2, flag]=sylvester3a(xx2,a,b,c,d);
        icount = icount + 1;
    end
    jcount = 0;
    for j=1:param_nbr,
        for i=j:param_nbr,
        jcount=jcount+1;
        x = xx2(:,:,jcount);
        elem1 = (D2g0(:,:,j,i)-D2g1(:,:,j,i)*A);
        elem1 = elem1 -( Dg1(:,:,j)*xx(:,:,i) + Dg1(:,:,i)*xx(:,:,j) );
        elemj0 = Dg0(:,:,j)-Dg1(:,:,j)*A-GAM1*xx(:,:,j);
        elemi0 = Dg0(:,:,i)-Dg1(:,:,i)*A-GAM1*xx(:,:,i);
        elem0 = elemj0*dB(:,:,i)+elemi0*dB(:,:,j);
        y = inva * (D2g3(:,:,j,i)-elem0-(elem1-GAM1*x)*B);
        d2B(:,:,j+length(indexo),i+length(indexo)) = y;
        d2B(:,:,i+length(indexo),j+length(indexo)) = y;
        y = y*M_.Sigma_e*B'+B*M_.Sigma_e*y'+ ...
            dB(:,:,j)*M_.Sigma_e*dB(:,:,i)'+dB(:,:,i)*M_.Sigma_e*dB(:,:,j)';
%         x = x(nauxe+1:end,nauxe+1:end);
%         y = y(nauxe+1:end,nauxe+1:end);
        d2A(:,:,j+length(indexo),i+length(indexo)) = x;
        d2A(:,:,i+length(indexo),j+length(indexo)) = x;
        d2Om(:,:,j+length(indexo),i+length(indexo)) = y;
        d2Om(:,:,i+length(indexo),j+length(indexo)) = y;
        end
    end    
    if ~isempty(indexo),
        d2Sig = zeros(M_.exo_nbr,M_.exo_nbr,length(indexo));
        for j=1:length(indexo)
            d2Sig(indexo(j),indexo(j),j) = 2;
            y = B*d2Sig(:,:,j)*B';
            d2Om(:,:,j,j) = y;
%             y = y(nauxe+1:end,nauxe+1:end);
%             H(:,j) = [zeros((m-nauxe)^2,1); dyn_vech(y)];
%             H(:,j) = [zeros(m^2,1); dyn_vech(y)];
%             dOm(:,:,j) = y;
            for i=1:param_nbr, 
                y = dB(:,:,i)*dSig(:,:,j)*B'+B*dSig(:,:,j)*dB(:,:,i)';
                d2Om(:,:,j,i+length(indexo)) = y;
                d2Om(:,:,i+length(indexo),j) = y;
            end
        end
    end
end

return

function g22 = get_2nd_deriv(gpp,m,n,i,j),

g22=zeros(m,n);
is=find(gpp(:,3)==i);
is=is(find(gpp(is,4)==j));

if ~isempty(is),
    g22(sub2ind([m,n],gpp(is,1),gpp(is,2)))=gpp(is,5)';
end
return

function g22 = get_all_2nd_derivs(gpp,m,n,npar),

g22=zeros(m,n,npar,npar);
% c=ones(npar,npar);
% c=triu(c);
% ic=find(c);

for is=1:length(gpp),
%     d=zeros(npar,npar);
%     d(gpp(is,3),gpp(is,4))=1;
%     indx = find(ic==find(d));
    g22(gpp(is,1),gpp(is,2),gpp(is,3),gpp(is,4))=gpp(is,5);
    g22(gpp(is,1),gpp(is,2),gpp(is,4),gpp(is,3))=gpp(is,5);
end

return

function r22 = get_all_resid_2nd_derivs(rpp,m,npar),

r22=zeros(m,npar,npar);
% c=ones(npar,npar);
% c=triu(c);
% ic=find(c);

for is=1:length(rpp),
%     d=zeros(npar,npar);
%     d(rpp(is,2),rpp(is,3))=1;
%     indx = find(ic==find(d));
    r22(rpp(is,1),rpp(is,2),rpp(is,3))=rpp(is,4);
    r22(rpp(is,1),rpp(is,3),rpp(is,2))=rpp(is,4);
end

return