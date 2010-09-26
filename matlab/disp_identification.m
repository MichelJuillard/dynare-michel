function disp_identification(pdraws, idemodel, idemoments, disp_pcorr)

% Copyright (C) 2008-2010 Dynare Team
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

global bayestopt_

if nargin<4 | isempty(disp_pcorr),
    disp_pcorr=0;
end

[SampleSize, npar] = size(pdraws);
jok = 0;
jokP = 0;
jokJ = 0;
jokPJ = 0;
if ~any(any(idemodel.ind==0))
    disp(['All parameters are identified in the model in the MC sample (rank of H).' ]),
    disp(' ')
end
if ~any(any(idemoments.ind==0))
    disp(['All parameters are identified by J moments in the MC sample (rank of J)' ]),
end
for j=1:npar,
    if any(idemodel.ind(j,:)==0),
        pno = 100*length(find(idemodel.ind(j,:)==0))/SampleSize;
        disp(['Parameter ',bayestopt_.name{j},' is not identified in the model for ',num2str(pno),'% of MC runs!' ])
        disp(' ')
    end
    if any(idemoments.ind(j,:)==0),
        pno = 100*length(find(idemoments.ind(j,:)==0))/SampleSize;
        disp(['Parameter ',bayestopt_.name{j},' is not identified by J moments for ',num2str(pno),'% of MC runs!' ])
        disp(' ')
    end
    if any(idemodel.ind(j,:)==1),
        iok = find(idemodel.ind(j,:)==1);
        jok = jok+1;
        kok(jok) = j;
        mmin(jok,1) = min(idemodel.Mco(j,iok));
        mmean(jok,1) = mean(idemodel.Mco(j,iok));
        mmax(jok,1) = max(idemodel.Mco(j,iok));
        [ipmax, jpmax] = find(abs(squeeze(idemodel.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
        if ~isempty(ipmax)
            jokP = jokP+1;
            kokP(jokP) = j;
            ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
            [N,X]=hist(ipmax,[1:npar]);
            jpM(jokP)={find(N)};
            NPM(jokP)={N(find(N))./SampleSize.*100};
            pmeanM(jokP)={mean(squeeze(idemodel.Pco(j,find(N),iok))')};
            pminM(jokP)={min(squeeze(idemodel.Pco(j,find(N),iok))')};
            pmaxM(jokP)={max(squeeze(idemodel.Pco(j,find(N),iok))')};
        end
    end
    if any(idemoments.ind(j,:)==1),
        iok = find(idemoments.ind(j,:)==1);
        jokJ = jokJ+1;
        kokJ(jokJ) = j;
        mminJ(jokJ,1) = min(idemoments.Mco(j,iok));
        mmeanJ(jokJ,1) = mean(idemoments.Mco(j,iok));
        mmaxJ(jokJ,1) = max(idemoments.Mco(j,iok));
        [ipmax, jpmax] = find(abs(squeeze(idemoments.Pco(j,[1:j-1,j+1:end],iok)))>0.95);
        if ~isempty(ipmax)
            jokPJ = jokPJ+1;
            kokPJ(jokPJ) = j;
            ipmax(find(ipmax>=j))=ipmax(find(ipmax>=j))+1;
            [N,X]=hist(ipmax,[1:npar]);
            jpJ(jokPJ)={find(N)};
            NPJ(jokPJ)={N(find(N))./SampleSize.*100};
            pmeanJ(jokPJ)={mean(squeeze(idemoments.Pco(j,find(N),iok))')};
            pminJ(jokPJ)={min(squeeze(idemoments.Pco(j,find(N),iok))')};
            pmaxJ(jokPJ)={max(squeeze(idemoments.Pco(j,find(N),iok))')};
        end
    end
end


dyntable('Multi collinearity in the model:',char('param','min','mean','max'), ...
         char(bayestopt_.name(kok)),[mmin, mmean, mmax],10,10,6);
disp(' ')
for j=1:npar,
    iweak = length(find(idemodel.Mco(j,:)'>(1-1.e-10)));
    if iweak, 
        disp('WARNING !!!')
        disp(['Model derivatives of parameter ',bayestopt_.name{j},' are multi-collinear (with tol = 1.e-10) for ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
        if npar>(j+1),
            [ipair, jpair] = find(squeeze(idemodel.Pco(j,j+1:end,:))'>(1-1.e-10));
        else
            [ipair, jpair] = find(squeeze(idemodel.Pco(j,j+1:end,:))>(1-1.e-10));
        end
        if ~isempty(jpair),
            for jx=j+1:npar,
                ixp = find(jx==(jpair+j));
                if ~isempty(ixp)
                    disp(['Model derivatives of parameters [',bayestopt_.name{j},',',bayestopt_.name{jx},'] are collinear (with tol = 1.e-10) for ',num2str(length(ixp)/SampleSize*100),'% of MC runs!' ])            
                end
            end
        end
    end
end
disp(' ')

dyntable('Multi collinearity for moments in J:',char('param','min','mean','max'), ...
         char(bayestopt_.name(kokJ)),[mminJ, mmeanJ, mmaxJ],10,10,6);
disp(' ')
for j=1:npar,
    iweak = length(find(idemoments.Mco(j,:)'>(1-1.e-10)));
    if iweak, 
        disp('WARNING !!!')
        disp(['Moment derivatives of parameter ',bayestopt_.name{j},' are multi-collinear (with tol = 1.e-10) for ',num2str(iweak/SampleSize*100),'% of MC runs!' ])
        if npar>(j+1),
            [ipair, jpair] = find(squeeze(idemoments.Pco(j,j+1:end,:))'>(1-1.e-10));
        else
            [ipair, jpair] = find(squeeze(idemoments.Pco(j,j+1:end,:))>(1-1.e-10));
        end
        if ~isempty(jpair),
            for jx=j+1:npar,
                ixp = find(jx==(jpair+j));
                if ~isempty(ixp)
                    disp(['Moment derivatives of parameters [',bayestopt_.name{j},',',bayestopt_.name{jx},'] are collinear (with tol = 1.e-10) for ',num2str(length(ixp)/SampleSize*100),'% of MC runs!' ])            
                end
            end
        end
    end
end
disp(' ')

if disp_pcorr,
    for j=1:length(kokP),
        dyntable([bayestopt_.name{kokP(j)},' pairwise correlations in the model'],char(' ','min','mean','max'), ...
                 char(bayestopt_.name{jpM{j}}),[pminM{j}' pmeanM{j}' pmaxM{j}'],10,10,3);  
    end

    for j=1:length(kokPJ),
        dyntable([bayestopt_.name{kokPJ(j)},' pairwise correlations in J moments'],char(' ','min','mean','max'), ...
                 char(bayestopt_.name{jpJ{j}}),[pminJ{j}' pmeanJ{j}' pmaxJ{j}'],10,10,3);  
    end
end
