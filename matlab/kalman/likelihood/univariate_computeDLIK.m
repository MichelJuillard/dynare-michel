function [Da,DP1,DLIK,D2a,D2P,Hesst] = univariate_computeDLIK(k,indx,Z,Zflag,v,K,PZ,F,Da,DYss,DP,DH,notsteady,D2a,D2Yss,D2P)

% Copyright (C) 2012 Dynare Team
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

% AUTHOR(S) marco.ratto@jrc.ec.europa.eu

persistent DDK DDF DD2K DD2F

if notsteady,
    if Zflag
        Dv   = -Z*Da(:,:) - Z*DYss(:,:);
        DF = zeros(k,1);
        DK = zeros([rows(K),k]);
        for j=1:k,
            DF(j)=Z*DP(:,:,j)*Z'+DH;
            DK(:,j) = (DP(:,:,j)*Z')/F-PZ*DF(j)/F^2;
        end
        if nargout>4
            D2F = zeros(k,k);
            D2v = zeros(k,k);
            D2K = zeros(rows(K),k,k);
            for j=1:k,
                D2v(:,j)   = -Z*D2a(:,:,j) - Z*D2Yss(:,:,j);
                for i=j:k,
                    D2F(j,i)=Z*D2P(:,:,j,i)*Z';
                    D2F(i,j)=D2F(j,i);
                    D2K(:,j,i) = (D2P(:,:,j,i)*Z')/F-(DP(:,:,j)*Z')*DF(i)/F^2-(DP(:,:,i)*Z')*DF(j)/F^2- ...
                        PZ*D2F(j,i)/F^2 + 2*PZ/F^3*DF(i)*DF(j);
                    D2K(:,i,j) = D2K(:,j,i);
                end
            end
        end
        
    else
        Dv   = -Da(Z,:) - DYss(Z,:);
        DF = squeeze(DP(Z,Z,:))+DH';
        DK = squeeze(DP(:,Z,:))/F-PZ*transpose(DF)/F^2;
        if nargout>4
            D2v   = squeeze(-D2a(Z,:,:) - D2Yss(Z,:,:));
            D2F = squeeze(D2P(Z,Z,:,:));
            D2K = squeeze(D2P(:,Z,:,:))/F;
            for j=1:k,
                D2K(:,:,j) = D2K(:,:,j) -PZ*D2F(j,:)/F^2 - squeeze(DP(:,Z,:))*DF(j)/F^2 - ...
                    squeeze(DP(:,Z,j))*DF'/F^2 + 2/F^3*PZ*DF'*DF(j);
            end
        end
    end
    if nargout>4
        DD2K(:,indx,:,:)=D2K;
        DD2F(indx,:,:)=D2F;
    end
    DDK(:,indx,:)=DK;
    DDF(indx,:)=DF;
else
    DK = squeeze(DDK(:,indx,:));
    DF = DDF(indx,:)';
    if nargout>4
        D2K = squeeze(DD2K(:,indx,:,:));
        D2F = squeeze(DD2F(indx,:,:));
    end
    if Zflag
        Dv   = -Z*Da(:,:) - Z*DYss(:,:);
        if nargout>4
            D2v = zeros(k,k);
            for j=1:k,
                D2v(:,j)   = -Z*D2a(:,:,j) - Z*D2Yss(:,:,j);
            end
        end
    else
        Dv   = -Da(Z,:) - DYss(Z,:);
        if nargout>4
            D2v   = squeeze(-D2a(Z,:,:) - D2Yss(Z,:,:));
        end
    end
end

DLIK = DF/F + 2*Dv'/F*v - v^2/F^2*DF;
if nargout==6
    Hesst = D2F/F-1/F^2*(DF*DF') + 2*D2v/F*v + 2*(Dv'*Dv)/F - 2*(DF*Dv)*v/F^2 ...
        - v^2/F^2*D2F - 2*v/F^2*(Dv'*DF') + 2*v^2/F^3*(DF*DF');
elseif nargout==4,
    D2a = 1/F^2*(DF*DF') + 2*(Dv'*Dv)/F ;
%     D2a = -1/F^2*(DF*DF') + 2*(Dv'*Dv)/F  + 2*v^2/F^3*(DF*DF') ...
%         - 2*(DF*Dv)*v/F^2 - 2*v/F^2*(Dv'*DF');
%     D2a = +2*(Dv'*Dv)/F + (DF' * DF)/F^2;
end

Da = Da + DK*v+K*Dv;
if nargout>4
    
    D2a = D2a + D2K*v;
    for j=1:k,
        %         D2a(:,:,j) = D2a(:,:,j) + DK*Dv(j) + DK(:,j)*Dv + K*D2v(j,:);
        for i=1:j,
            D2a(:,j,i) = D2a(:,j,i) + DK(:,i)*Dv(j) + DK(:,j)*Dv(i) + K*D2v(j,i);
            D2a(:,i,j) = D2a(:,j,i);
        end
    end
    
    
end

if notsteady,
    DP1 = DP*0;
    if Zflag,
        for j=1:k,
            DP1(:,:,j)=DP(:,:,j) - (DP(:,:,j)*Z')*K'-PZ*DK(:,j)';
        end
    else
        for j=1:k,
            DP1(:,:,j)=DP(:,:,j) - (DP(:,Z,j))*K'-PZ*DK(:,j)';
        end
    end
    if nargout>4,
        if Zflag,
            for j=1:k,
                D2P = D2P;
            end
        else
            D2PZ = squeeze(D2P(:,Z,:,:));
            DPZ = squeeze(DP(:,Z,:));
            for j=1:k,
                for i=j:k,
                    D2P(:,:,j,i) = D2P(:,:,j,i) - D2PZ(:,j,i)*K' - DPZ(:,j)*DK(:,i)'- DPZ(:,i)*DK(:,j)' - PZ*squeeze(D2K(:,j,i))';
                    D2P(:,:,i,j) = D2P(:,:,j,i);
                end
            end
            
        end
    end
else
    DP1=DP;
end


