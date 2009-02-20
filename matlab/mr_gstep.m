function [f0, x, ig] = mr_gstep(func0,x,htol0,varargin)
% function [f0, x] = mr_gstep(func0,x,htol0,varargin)
% 
% Gibbs type step in optimisation

% Copyright (C) 2006 Dynare Team
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

global bayestopt_ options_
persistent h1 

gstep_ = options_.gstep;
if nargin<3, 
    htol = 1.e-6;
else
    htol = htol0;
end
func = str2func(func0);
f0=feval(func,x,varargin{:});
n=size(x,1);
h2=bayestopt_.ub-bayestopt_.lb;

if isempty(h1),
    h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/4);
end

xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;
%for i=1:n,
i=0;
ig=zeros(n,1);
while i<n,
    i=i+1;
    h10=h1(i);
    hcheck=0;
    dx=[];
    xh1(i)=x(i)+h1(i);
    fx = feval(func,xh1,varargin{:});
    it=1;
    dx=(fx-f0);
    ic=0;
%     if abs(dx)>(2*htol),
%         c=mr_nlincon(xh1,varargin{:});
%         while c
%             h1(i)=h1(i)*0.9;
%             xh1(i)=x(i)+h1(i);
%             c=mr_nlincon(xh1,varargin{:});        
%             ic=1;
%         end   
%         if ic,
%             fx = feval(func,xh1,varargin{:});
%             dx=(fx-f0);
%         end
%     end
    
    icount = 0;
    h0=h1(i);
    while (abs(dx(it))<0.5*htol | abs(dx(it))>(2*htol)) & icount<10 & ic==0,
        %while abs(dx(it))<0.5*htol & icount< 10 & ic==0,
        icount=icount+1;
        if abs(dx(it)) ~= 0,
            if abs(dx(it))<0.5*htol
                h1(i)=min(0.3*abs(x(i)), 0.9*htol/abs(dx(it))*h1(i));
                xh1(i)=x(i)+h1(i);
%                 c=mr_nlincon(xh1,varargin{:});
%                 while c
%                     h1(i)=h1(i)*0.9;
%                     xh1(i)=x(i)+h1(i);
%                     c=mr_nlincon(xh1,varargin{:});        
%                     ic=1;
%                 end  
            end
            if abs(dx(it))>(2*htol),
                h1(i)= htol/abs(dx(it))*h1(i);
                xh1(i)=x(i)+h1(i);
            end
            fx = feval(func,xh1,varargin{:});
            it=it+1;
            dx(it)=(fx-f0);
            h0(it)=h1(i);
            if h1(i)<1.e-12*min(1,h2(i)),
                ic=1;
                hcheck=1;
            end
        else
            h1(i)=1;
            ic=1;
        end
    end
    f1(:,i)=fx;
    xh1(i)=x(i)-h1(i);
%     c=mr_nlincon(xh1,varargin{:});
%    ic=0;
%     while c
%         h1(i)=h1(i)*0.9;
%         xh1(i)=x(i)-h1(i);
%         c=mr_nlincon(xh1,varargin{:});  
%         ic = 1;
%     end    
    fx = feval(func,xh1,varargin{:});
    f_1(:,i)=fx;
%     if ic,
%         xh1(i)=x(i)+h1(i);
%         f1(:,i)=feval(func,xh1,varargin{:});
%     end
    if hcheck & htol<1,
        htol=min(1,max(min(abs(dx))*2,htol*10));
        h1(i)=h10;
        xh1(i)=x(i);
        i=i-1;
    else
        gg=zeros(size(x));    
        hh=gg;
        gg(i)=(f1(i)'-f_1(i)')./(2.*h1(i));
        if abs(f1(i)+f_1(i)-2*f0)>1.e-12,
            hh(i) = abs(1/( (f1(i)+f_1(i)-2*f0)./(h1(i)*h1(i)) ));
        else
            hh(i) = 1;
        end
            
        if gg(i)*(hh(i)*gg(i))/2 > htol,
            [f0 x fc retcode] = csminit(func0,x,f0,gg,0,diag(hh),varargin{:});
            ig(i)=1;
        end
        xh1=x;
    end
    save gstep
end

save gstep



