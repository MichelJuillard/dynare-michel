function [f0, x, ig] = mr_gstep(h1,x,func0,htol0,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% function [f0, x, ig] = mr_gstep(h1,x,func0,htol0,varargin)
%
% Gibbs type step in optimisation

% Copyright (C) 2006-2012 Dynare Team
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

n=size(x,1);
if isempty(h1),
    h1=DynareOptions.gradient_epsilon*ones(n,1);
end


if isempty(htol0)
    htol = 1.e-6;
else
    htol = htol0;
end
f0=feval(func0,x,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults);

xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;

i=0;
ig=zeros(n,1);
while i<n
    i=i+1;
    h10=h1(i);
    hcheck=0;
    dx=[];
    xh1(i)=x(i)+h1(i);
    fx = feval(func0,xh1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults);
    f1(:,i)=fx;
    xh1(i)=x(i)-h1(i);
    fx = feval(func0,xh1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults);
    f_1(:,i)=fx;
    if hcheck && htol<1
        htol=min(1,max(min(abs(dx))*2,htol*10));
        h1(i)=h10;
        xh1(i)=x(i);
        i=i-1;
    else
        gg=zeros(size(x));
        hh=gg;
        gg(i)=(f1(i)'-f_1(i)')./(2.*h1(i));
        hh(i) = 1/max(1.e-9,abs( (f1(i)+f_1(i)-2*f0)./(h1(i)*h1(i)) ));
        if gg(i)*(hh(i)*gg(i))/2 > htol
            [f0 x fc retcode] = csminit1(func0,x,f0,gg,0,diag(hh),DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults);
            ig(i)=1;
            fprintf(['Done for param %s = %8.4f\n'],BayesInfo.name{i},x(i))
        end
        xh1=x;
    end
    save gstep.mat x h1 f0
end

save gstep.mat x h1 f0

