function [x,%check]=solve(func,x,varargin)
%check=[];
// Copyright (C) 2001 Michel Juillard
// 
 

nn = size(x,1);

fjac = zeros(nn,nn);
g = zeros(nn,1);
 
tolf = %eps^(2/3);
tolmin = .000000000037;
tolx = .000000000037;
 
stpmx = 100;
maxit = 2000;
 
%check = 0;

if argn(2) > 2 
 fvec = evstr(func+'(x,varargin)');
else
  fvec = evstr(func+'(x)');
end

f = .5*fvec'*fvec;
 
 
if max(abs(fvec))<.01*tolf then
  return
end
 
stpmax = stpmx*max([abs(sqrt(x'*x));nn]);
if argn(2) > 2
 func_call = func+'(xdh,varargin)';
else
 func_call = func+'(xdh)';
end
 
for its = 1:maxit
  dh = max(abs(x),gstep_*ones(nn,1))*(%eps^(1/3));
  for j = 1:nn
    xdh = x;
    xdh(j) = xdh(j)+dh(j);
    fjac(:,j) = (evstr(func_call)-fvec) ./ dh(j);
    g(j) = fvec'*fjac(:,j);
  end
   
  if debug_ then
    cond(fjac);
  end
   
  p = -fjac\fvec;
  xold = x;
  fold = f;
   
  if argn(2) > 2
    [x,f,fvec,%check] = lnsrch(xold,fold,g,p,stpmax,func,varargin);
  else
    [x,f,fvec,%check] = lnsrch(xold,fold,g,p,stpmax,func);
  end 

  if %check>0 then
    den = max([f;.5*nn]);
    if max(abs(g) .* max([abs(x');ones(1,nn)])')/den<tolmin then
      return
       
    else
      dyn_disp(' ');
      dyn_disp('SOLVE: Iteration '+string(its));
      dyn_disp('Spurious convergence.');
      dyn_disp(x);
      return
       
    end
     
    if max(abs(x-xold) ./ max([abs(x);ones(1,nn)])')<tolx then
      dyn_disp(' ');
      dyn_disp('SOLVE: Iteration '+string(its));
      dyn_disp('Convergence on dX.');
      dyn_disp(x);
      return
       
    end
  elseif max(abs(fvec)) < tolf then
    return
     
  end
end
 
%check = 1;
dyn_disp(' ');
dyn_disp('SOLVE: maxit has been reached');
 
// 01/14/01 MJ lnsearch is now a separate function
// 01/16/01 MJ added varargin to function evaluation
// 04/13/01 MJ added test  f < tolf !!
// 05/11/01 MJ changed tests for 'check' so as to remove 'continue' which is
//             an instruction which appears only in version 6
// 09/26/01 MJ translated to Scilab
// 10/08/01 MJ changed end criterium from f < tolf to max(abs(fvec)) < tolf 
 
 
 
 
 
 
 
 
