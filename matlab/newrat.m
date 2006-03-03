function [xparam1, hh, gg, fval, igg] = newrat(func0, x, hh, gg, igg, ftol0, nit, flagg, varargin)
%
%  Copyright (C) 2004 Marco Ratto
%
%  [xparam1, hh, gg, fval, igg] = newrat(func0, x, hh, gg, igg, ftol0, nit, flagg, varargin)
%
%  Optimiser with outer product gradient and 'Gibbs type' steps
%  uses Chris Sims subroutine for line search
%
%  func0 = name of the function
%  there must be a version of the function called [func0,'_hh.m'], that also
%  gives as second OUTPUT the single contributions at times t=1,...,T
%    of the log-likelihood to compute outer product gradient
%
%  x = starting guess
%  hh = initial Hessian [OPTIONAL]
%  gg = initial gradient [OPTIONAL]
%  igg = initial inverse Hessian [OPTIONAL]
%  ftol0 = ending criterion for function change 
%  nit = maximum number of iterations
%
%  In each iteration, Hessian is computed with outer product gradient.
%  for final Hessian (to start Metropolis):
%  flagg = 0, final Hessian computed with outer product gradient
%  flagg = 1, final 'mixed' Hessian: diagonal elements computed with numerical second order derivatives
%             with correlation structure as from outer product gradient, 
%  flagg = 2, full numerical Hessian
%
%  varargin = list of parameters for func0 
%

  global bayestopt_
icount=0;
nx=length(x);
xparam1=x;
%ftol0=1.e-6;
flagit=0;  % mode of computation of hessian in each iteration
ftol=ftol0;
gtol=1.e-3;
htol=ftol0;
htol0=ftol0;

func_hh = [func0,'_hh'];
func = str2func(func0);
fval0=feval(func,x,varargin{:});
fval=fval0;
if isempty(hh)
    [dum, gg, htol0, igg, hhg]=mr_hessian(func_hh,x,flagit,htol,varargin{:});
    hh0 = reshape(dum,nx,nx);
    hh=hhg;
    if min(eig(hh0))<0,
        hh0=hhg; %generalized_cholesky(hh0);
    elseif flagit==2,
        hh=hh0;
        igg=inv(hh);
    end
    if htol0>htol,
        htol=htol0;
        ftol=htol0;
    end
else
    hh0=hh;
    hhg=hh;
    igg=inv(hh);
end
disp(['Gradient norm ',num2str(norm(gg))])
ee=eig(hh);
disp(['Minimum Hessian eigenvalue ',num2str(min(ee))])
disp(['Maximum Hessian eigenvalue ',num2str(max(ee))])
g=gg;
check=0;
if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
save m1 x hh g hhg igg fval0

igrad=1;
igibbs=1;
inx=eye(nx);
jit=0;
while norm(gg)>gtol & check==0 & jit<nit,
    jit=jit+1;
    tic
    icount=icount+1;
    bayestopt_.penalty = fval0(icount);
    disp([' '])
    disp(['Iteration ',num2str(icount)])
    [fval x0 fc retcode] = csminit(func0,xparam1,fval0(icount),gg,0,igg,varargin{:});
    if igrad,
        [fval1 x01 fc retcode1] = csminit(func0,xparam1,fval0(icount),gg,0,inx,varargin{:});
        if fval1<fval,
            fval=fval1;
            x0=x01;        
            disp('Gradient step!!')
        else
            igrad=0;
        end
    end
    if (fval0(icount)-fval)<1.e-2*(gg'*(igg*gg))/2 & igibbs,
        [fvala, x0] = mr_gstep(func0,x0,htol,varargin{:});
         if (fval-fvala)<5*(fval0(icount)-fval),
             igibbs=0;
             disp('Last Gibbs step, gain too small!!')
         else
            disp('Gibbs step!!')
        end
        fval=fvala;
    end
    if (fval0(icount)-fval)<ftol & flagit==0,
        disp('Try diagonal Hessian')
        ihh=diag(1./(diag(hhg)));        
        [fval2 x02 fc retcode2] = csminit(func2str(func),xparam1,fval0(icount),gg,0,ihh,varargin{:});
        if fval2<fval,
            x0=x02;
            fval=fval2;
            if (fval0(icount)-fval2)>=ftol ,
                %hh=diag(diag(hh));
                disp('Diagonal Hessian successful')            
            end
        end
    end        
    if (fval0(icount)-fval)<ftol & flagit==0,
        disp('Try gradient direction')
        ihh0=inx.*1.e-4;        
        [fval3 x03 fc retcode3] = csminit(func2str(func),xparam1,fval0(icount),gg,0,ihh0,varargin{:});
        if fval3<fval,
            x0=x03;
            fval=fval3;
            if (fval0(icount)-fval3)>=ftol ,
                %hh=hh0;
                %ihh=ihh0;
                disp('Gradient direction successful')            
            end
        end
    end        
    xparam1=x0;
    x(:,icount+1)=xparam1;
    fval0(icount+1)=fval;
    %if (fval0(icount)-fval)<ftol*ftol & flagg==1;,
    if (fval0(icount)-fval)<ftol,
        disp('No further improvement is possible!')
        check=1;
        if flagit==2,
            hh=hh0;
        elseif flagg>0,
            [dum, gg, htol0, igg, hhg]=mr_hessian(func_hh,xparam1,flagg,ftol0,varargin{:});   
            if flagg==2,
                hh = reshape(dum,nx,nx);
                ee=eig(hh);
                if min(ee)<0
                    hh=hhg;
                end
            else
                hh=hhg;
            end
        end
        disp(['Actual dxnorm ',num2str(norm(x(:,end)-x(:,end-1)))])
        disp(['FVAL          ',num2str(fval)])
        disp(['Improvement   ',num2str(fval0(icount)-fval)])
        disp(['Ftol          ',num2str(ftol)])
        disp(['Htol          ',num2str(htol0)])
        disp(['Gradient norm  ',num2str(norm(gg))])
        ee=eig(hh);
        disp(['Minimum Hessian eigenvalue ',num2str(min(ee))])
        disp(['Maximum Hessian eigenvalue ',num2str(max(ee))])
         g(:,icount+1)=gg;
    else
        
        df = fval0(icount)-fval;
        disp(['Actual dxnorm ',num2str(norm(x(:,end)-x(:,end-1)))])
        disp(['FVAL          ',num2str(fval)])
        disp(['Improvement   ',num2str(df)])
        disp(['Ftol          ',num2str(ftol)])
        disp(['Htol          ',num2str(htol0)])

        if df<htol0,
            htol=max(ftol0,df/10);
        end
        
        if norm(x(:,icount)-xparam1)>1.e-12,
            save m1 x fval0 -append
            [dum, gg, htol0, igg, hhg]=mr_hessian(func_hh,xparam1,flagit,htol,varargin{:});
            if htol0>ftol,
                ftol=htol0;
                htol=htol0;
                disp(' ')
                disp('Numerical noise in the likelihood')
                disp('Tolerance has to be relaxed')
                disp(' ')
            elseif htol0<ftol,
                ftol=max(htol0, ftol0);
            end
            hh0 = reshape(dum,nx,nx);
            hh=hhg;
            if flagit==2,
                if min(eig(hh0))<=0,
                    hh0=hhg; %generalized_cholesky(hh0);
                else 
                    hh=hh0;
                    igg=inv(hh);
                end
            end
        end
        disp(['Gradient norm  ',num2str(norm(gg))])
        ee=eig(hh);
        disp(['Minimum Hessian eigenvalue ',num2str(min(ee))])
        disp(['Maximum Hessian eigenvalue ',num2str(max(ee))])
        if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
        t=toc;
        disp(['Elapsed time for iteration ',num2str(t),' s.'])
        
         g(:,icount+1)=gg;
        save m1 x hh g hhg igg fval0
    end
end

save m1 x hh g hhg igg fval0
if ftol>ftol0,
    disp(' ')
    disp('Numerical noise in the likelihood')
    disp('Tolerance had to be relaxed')
    disp(' ')
end

if jit==nit,
    disp(' ')
    disp('Maximum number of iterations reached')
    disp(' ')
end

if norm(gg)<=gtol,
    disp(['Estimation ended:'])
    disp(['Gradient norm < ', num2str(gtol)])
end
if check==1,
    disp(['Estimation successful.'])
end

return

%  
function f00 = lsearch(lam,func,x,dx,varargin)


x0=x-dx*lam;
f00=feval(func,x0,varargin{:});






