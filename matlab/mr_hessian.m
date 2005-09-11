% Copyright (C) 2001 Michel Juillard
%
% computes second order partial derivatives
% uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
%
% Adapted by M. Ratto from original M. Juillard routine
%
function [hessian_mat, gg] = mr_hessian(func,x,hflag,varargin)
global gstep_
persistent h1

func = str2func(func);
f0=feval(func,x,varargin{:});
n=size(x,1);
%h1=max(abs(x),gstep_*ones(n,1))*eps^(1/3);
%h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/6);
if isempty(h1),
    h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/4);
end
htol = 1.e-4;
xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;

for i=1:n,
    xh1(i)=x(i)+h1(i);
    fx=feval(func,xh1,varargin{:});
    if abs(fx-f0)<htol | abs(fx-f0)>(5*htol),
        c=mr_nlincon(xh1,varargin{:});
        ic=0;
        while c
            h1(i)=h1(i)*0.9;
            xh1(i)=x(i)+h1(i);
            c=mr_nlincon(xh1,varargin{:});        
            ic=1;
        end   
        if ic,
            fx=feval(func,xh1,varargin{:});
        end

        icount = 0;
        while abs(fx-f0)<htol & icount< 10,
            icount=icount+1;
            h1(i)=min(0.3*abs(x(i)), 1.e-3/(abs(fx-f0)/h1(i)));
            xh1(i)=x(i)+h1(i);
            c=mr_nlincon(xh1,varargin{:});
            while c
                h1(i)=h1(i)*0.9;
                xh1(i)=x(i)+h1(i);
                c=mr_nlincon(xh1,varargin{:});        
            end        
            fx=feval(func,xh1,varargin{:});
        end
        while abs(fx-f0)>(5*htol),
            h1(i)=h1(i)*0.5;    
            xh1(i)=x(i)+h1(i);
            fx=feval(func,xh1,varargin{:});
        end
    end
    f1(:,i)=fx;
    xh1(i)=x(i)-h1(i);
    c=mr_nlincon(xh1,varargin{:});
    ic=0;
    while c
        h1(i)=h1(i)*0.9;
        xh1(i)=x(i)-h1(i);
        c=mr_nlincon(xh1,varargin{:});  
        ic = 1;
    end    
    fx=feval(func,xh1,varargin{:});
    f_1(:,i)=fx;
    if ic,
        xh1(i)=x(i)+h1(i);
        f1(:,i)=feval(func,xh1,varargin{:});
    end
    xh1(i)=x(i);
end

h_1=h1;
xh1=x;
xh_1=xh1;
gg=(f1'-f_1')./(2.*h1);

if hflag,
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n
        if i > 1
            k=[i:n:n*(i-1)];
            hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
        end 
        hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        temp=f1+f_1-f0*ones(1,n);
        for j=i+1:n
            xh1(i)=x(i)+h1(i);
            xh1(j)=x(j)+h_1(j);
            xh_1(i)=x(i)-h1(i);
            xh_1(j)=x(j)-h_1(j);
            %hessian_mat(:,(i-1)*n+j)=-(-feval(func,xh1,varargin{:})-feval(func,xh_1,varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
            temp1 = feval(func,xh1,varargin{:});
            %c=mr_nlincon(xh1,varargin{:});
            %if c, disp( ['hessian warning cross ', num2str(c) ]), end
            
            temp2 = feval(func,xh_1,varargin{:});
            %c=mr_nlincon(xh_1,varargin{:});
            %if c, disp( ['hessian warning cross ', num2str(c) ]), end
            hessian_mat(:,(i-1)*n+j)=-(-temp1 -temp2+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
            xh1(i)=x(i);
            xh1(j)=x(j);
            xh_1(i)=x(i);
            xh_1(j)=x(j);
            j=j+1;
        end
        i=i+1;
    end
    
else
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n,
        dum = (f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        if dum>0,
            hessian_mat(:,(i-1)*n+i)=dum;
        else
            hessian_mat(:,(i-1)*n+i)=gg(i)^2;
        end                        
    end
end
hh1=h1;
save hess
% 11/25/03 SA Created from Hessian_sparse (removed sparse)


