function [Xkk,Pkk,er,ykk1,Xkk1,yp,ep,ykk,XkN,ys,PkN,ers]=...
    ksirw2c(X,RWtype,Miss,Interv,F,P0,Q,iout,x0,nf,intD,QPflag)
% [Xkk,Pkk,er,ykk1,Xkk1,yp,ep,ykk,XkN,ys,PkN,ers]=...
%                ksirw2c(X,RWtype,Miss,Interv,F,P0,Q,iout,x0,nf,intD,QPflag)
%                        1   2      3    4    5 6  7  8   9  10  11   12
% Fast version of KALMSMO for use with SDR  
% No interventions, no missing values, only RW and IRW single parameter models
% Execution time gain - up to 10 times.
%  WT, September 2004

% Input:
%   X: (Ndat,Np+1)   - data (y-s in the first column)
%   RWtype (Np) - number of aux. parameters (eg. slopes)
%                 for given parameters
%   Miss: (Ndat,1)   - missing data flags (1 - missing sample)
%   Interv: (Ndat,1) - intervention flags (1 - intervention)
%   F  : (sqare matrix Np2) - transition matrix
%   Q  : (sqare matrix Np2) - noise variance (nvr) matrix
%   P0 : (sqare matrix Np2) - parameters cov. matrix (Pkk) - initial
%   iout    (integer)       - intermediate results on/off
%   x0                      - initial cond. for KF 
%   nf                      - max forecasting horizon 
%   intD                    - diagonal of Q matrix defining the 
%                             intervention type (size of diag(Q))
% Output:
%    XkN: as Xkk but smoothed
%    ys  : smoothed output 
%    PkN: (square matrix Np2 repeated 1:Ndat) - smoothed par. cov. mtx
%    Xkk: (Ndat,Np2)         - filtered parameters history
%    Pkk :                   - P(k|k) history
%    er : (Ndat)             - normalised error
%    ykk1: (Ndat)            - one step ahead prediction
%    yp: (Ndat,nf-1)         - predictions: yp(k,j)=y(k+j|k-1)
%    ep: (Ndat,nf-1)         - nrm. fact.of prediction errors: 
%                            -                       ep(k|j)=s(k+j|k-1)
%    ers : (Ndat)            - normalised error (smoothing)
% Q and P algorithms implemented (QPflag==1 : Q alg)

% Copyright (c) 2004 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Revision history :
%   20/09/2004, WT, converted to ksrivc from kalmsmo

% simplifications:
% RWtype = 1 or 0
% F = [1 1;0 1] or 1
% Q = [0 0;0 nvr] or nvr
% no interventions
% no missing data
% no forecasting built in

yp=[];ep=[];ykk=[];XkN=[];ys=[];PkN=[];ers=[];

% sizes and initial conds.
[Ndat Np]=size(X);Np=Np-1;Np2=length(F);

P=P0;

% defaults
if nargin<12
    QPflag=1;
end


if (nargin<9); x0=zeros(Np2,1); end
if isempty(x0),x0=zeros(Np2,1); end


if nargin<8, iout=0; end  % iout
if isempty(iout), iout=1; end

er=zeros(Ndat,1);
perr=zeros(Ndat,1);
Xkk=zeros(Np2,Ndat); Xkk1=[]; % Xkk1=Xkk; 

Pkk=zeros(Np2,Np2*Ndat); Pkk(:,1:Np2)=P0;

ykk1=zeros(Ndat,1); ykk=ykk1;
PHmat=zeros(Np2,Ndat);

nvr=Q(Np2,Np2);
Xhat=zeros(Np2,1);

if Np2==2   %  ****************************  IRW
    
    if ~isempty(x0), Xhat(1)=x0(1);Xhat(2)=x0(2);end
    
    for it=1:Ndat 
        u=X(it,2);
        % KF prediction step
        Xhat=[ Xhat(1)+Xhat(2); Xhat(2)];    % state X(k|k-1)
        P=[ P(1,1)+P(2,1)+P(1,2)+P(2,2) P(1,2)+P(2,2);P(2,1)+P(2,2) P(2,2)+nvr];  % cov. matrix
        
        HP=[ u*P(1,1), u*P(1,2)];
        div=1+u*P(1,1)*u;
        er(it)=div;div=1/div;
        yhat=u*Xhat(1);ykk1(it)=yhat;          % calc. and save save y(k|k-1)
        perr(it)=X(it,1)-yhat;                 % innovation

        % KF correction step
        Xhat=Xhat+HP'*perr(it)*div;            % X(k|k)
        yhat=u*Xhat(1);                        % yhat(k|k)
        P=P-u*u*[P(1,1)*P(1,1) P(1,1)*P(1,2);P(1,2)*P(1,1) P(1,2)*P(1,2)]*div;
        ykk(it)=yhat;                        % save 
        Xkk(:,it)=Xhat; 
        Pkk(:,(it-1)*Np2+1:it*Np2)=P;
        PHmat(:,it)=u*[ P(1,1);P(2,1)];
    end;                                   % of filtering loop
    
    
    L=zeros(2,1); 
    XkN=zeros(Np2,Ndat);   XkN(:,Ndat)=Xhat;
    ys=zeros(Ndat,1);      ys(Ndat)=Xhat(1);
    PkN=zeros(Np2,Np2*Ndat);
    PkN(:,(Ndat-1)*Np2+1:Ndat*Np2)=P;
    ers=zeros(Ndat,1);
    sP=P; 
    ers(Ndat)=er(Ndat);
    
    % smoothing loop
    for it=Ndat:-1:1               
        u=X(it,2);
        XkN(:,it)=Xhat;
        if QPflag, 
            ys(it)=u*Xhat(1);
        else
            pkk=Pkk(:,(it-1)*Np2+1:(it)*Np2);
            Xhat=Xkk(:,it)-[(pkk(1,1)+pkk(1,2))*L(1)+pkk(1,2)*L(2); (pkk(2,1)+pkk(2,2))*L(1)+pkk(2,2)*L(2)];
        end
        
%        if Miss(it)==0
            L=[(1-PHmat(1,it)*u)*(L(1)-perr(it)*u)-PHmat(2,it)*u*(L(1)+L(2)); L(1)+L(2)];
%        else
%            L=[L(1);L(1)+L(2)];
%        end
        if QPflag,
            Xhat=[Xhat(1)-Xhat(2)-nvr*L(2);Xhat(2)+nvr*L(2)];
        else
            ys(it)=u*Xhat(1);
            XkN(:,it)=Xhat;
        end
        if (it<Ndat)&(nargout>10)    % calculation of PkN requested 
                P=Pkk(:,(it-1)*Np2+1:it*Np2);
                Pdd=(P(1,1)*P(2,2)+P(1,1)*nvr+P(2,1)*nvr+P(1,2)*nvr+P(2,2)*nvr-P(2,1)*P(1,2));
                Pcc=(-P(1,1)*P(2,2)+P(2,1)*P(1,2));
                Paa=(Pcc/Pdd);
                Pbb=(sP(1,1)-P(1,1)-P(2,1)-P(1,2)-P(2,2));
                Pff=(sP(2,1)-P(2,1)-P(2,2));
                Pgg=(nvr*(P(2,1)+P(2,2))/Pdd*(sP(1,2)-P(1,2)-P(2,2))-Paa*(sP(2,2)-P(2,2)-nvr));
                Phh=(P(1,1)*P(2,2)+P(1,1)*nvr+P(1,2)*nvr-P(2,1)*P(1,2));
                Pbd=Pbb/Pdd;
                Phd=(Phh/Pdd);
                Pah=(sP(2,2)-P(2,2)-nvr);
                Pah1=Phd*(sP(1,2)-P(1,2)-P(2,2));
                Ph0=(Phh*Pbd+Paa*Pff);
                Pn=nvr*(P(2,1)+P(2,2));
                sP=[P(1,1)+Ph0*Phd+(Pah1+Paa*Pah)*Paa     P(1,2)+Ph0*(Pn/Pdd)-(Pah1+Paa*Pah)*Paa;
                    P(2,1)+(Pn*Pbd-Paa*Pff)*Phd+Pgg*Paa   P(2,2)+(Pn*Pbd-Paa*Pff)*(Pn/Pdd)-Pgg*Paa];
                PkN(:,(it-1)*Np2+1:it*Np2)=sP;
            ers(it)=1+u*u*sP(1,1);
        end
        
        if(abs(iout)>2)&~rem(it,10)            % display current est. 
            disp(Xhat(ipp)'); sPd=diag(sP)'; 
            disp(sqrt(sPd(ipp))) 
        end
    end; 
    
    ys=ys';   % consistent with previous versions
    
else % if Np2==2    ************************************   RW
    
    if ~isempty(x0), Xhat=x0;end   
    for it=1:Ndat 
        u=X(it,2);
        Xhat=Xhat;
        P=P+nvr;
        HP=u*P;
        div=1+u*P*u;
        er(it)=div;div=1/div;
        yhat=u*Xhat;ykk1(it)=yhat;          % calc. and save save y(k|k-1)
        perr(it)=X(it,1)-yhat;
        Xhat=Xhat+HP*perr(it)*div;
        yhat=u*Xhat;
        P=P-u*u*P*P*div;
        ykk(it)=yhat;                        % save 
        Xkk(:,it)=Xhat; 
        Pkk(it)=P;
        PHmat(it)=u*P;
    end;                                   % of filtering loop
    
    L=0; 
    XkN=zeros(Np2,Ndat);   XkN(:,Ndat)=Xhat;
    ys=zeros(Ndat,1);      ys(Ndat)=Xhat(1);
    PkN=zeros(Np2,Np2*Ndat);
    PkN(:,(Ndat-1)*Np2+1:Ndat*Np2)=P;
    ers=zeros(Ndat,1);
    sP=P; 
    ers(Ndat)=er(Ndat);
    for it=Ndat:-1:1               
        u=X(it,2);
        XkN(:,it)=Xhat;
        if QPflag, 
            ys(it)=u*Xhat;
        else
            pkk=Pkk(:,(it-1)*Np2+1:(it)*Np2);
            Xhat=Xkk(it)-Pkk(it)*L;
        end
        
%        if Miss(it)==0
            L=(1-u*PHmat(it))*(L-u*perr(it));
%        else
%            L=L;                           
%        end
        if QPflag,
            Xhat=Xhat+nvr*L;
        else
            ys(it)=u*Xhat;
            XkN(it)=Xhat;
        end
        if (it<Ndat)&(nargout>10)     
            P=Pkk(it);
            Pp=P+Q;
            PT=P/Pp;
            sP=P+PT*(sP-Pp)*PT;
            PkN(it)=sP;
            ers(it)=1+sP;
        end
    end; 
    
    ys=ys';   % consistent with previous versions    
    
end  % if Np2==2    
% end of m-file