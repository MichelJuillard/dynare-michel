function [J]=diffext(f,x,options,varargin)
[nargout,nargin] = argn(0)
//---------------------------------------------------------------------------
//DIFFEXT   Numerical approximation for hessian.
//          The method is Richardson`s extrapolation.
// Sample call
//   [D,err,relerr,n] = diffext('f',x,delta,toler)
// Inputs
//   f        name of the function
//   x        differentiation point
//   options  matrix of algorithm parameters
//   delta    error goal (1e-12) (suppressed MJ 02/27/02)
//   toler    relative error goal (1e-12)
// Return
//   J        Jacobian
// 
// NUMERICAL METHODS: MATLAB Programs, (c) John H. Mathews 1995
// 
// Modified F. Collard, August 2001
//---------------------------------------------------------------------------
if nargin > 2 then
  if ~(options==[]) then
//    %delta = options(1);
    toler = options(2);
    maxit = options(3);
  else
//    %delta = 1e-12;
    toler = 1e-12;
    maxit = 20;
  end
else
  %delta = 1e-12;
  toler = 1e-12;
  maxit = 20;
end
con = 1.4;
con2 = con*con;
big = 1e30;
safe = 2;
if nargin > 3
  ff = evstr(f+'(x,varargin)');
else
  ff = evstr(f+'(x)');
end
nx = size(x,1);
nf = size(ff,1);
J = zeros(nf,nx);
for xi = 1:nx
  err = big*ones(nf,1);;
//  relerr = big*ones(nf,1);
  h = max(abs(x(xi))/10, 10*gstep_)*100*gstep_;
  dx = zeros(nx,1);
  dx(xi,1) = h;
  if nargin > 3
    fs = evstr(f+'(x+dx,varargin)');
    fm = evstr(f+'(x-dx,varargin)');
  else
    fs = evstr(f+'(x+dx)');
    fm = evstr(f+'(x-dx)');
  end      
  D1 = (fs-fm)/(2*h);
  mask = ones(nf,1);
  j = 2;
  while (1)
    D2 = zeros(nf,j);
    h = h/con;
    dx(xi,1) = h;
    if nargin > 3
      fs = evstr(f+'(x+dx,varargin)');
      fm = evstr(f+'(x-dx,varargin)');
    else
      fs = evstr(f+'(x+dx)');
      fm = evstr(f+'(x-dx)');
    end      
    for fi = 1:nf
      if mask(fi)
	err2 = big;
	D2(fi,1) = (fs(fi)-fm(fi))/(2*h);
	fac = con2;
	for k = 2:j
	  D2(fi,k) = D2(fi,k-1)+(D2(fi,k-1)-D1(fi,k-1))/(fac-1);
	  fac = con2*fac;
	  errt = max(abs(D2(fi,k)-D2(fi,k-1)),abs(D2(fi,k)-D1(fi,k-1)));
	  if errt <= err2
	    err2 = errt;
	    deriv = D2(fi,k); 
	  end
	end
	err(fi) = abs(D2(fi,j)-D1(fi,j-1));
//	relerr(fi) = 2*err(fi)/(abs(D2(fi,j))+abs(D1(fi,j-1))+%eps);
//	if (err(fi)  < toler & relerr(fi) < %delta)| err(fi) > safe*err2
	if err(fi)  < toler | err(fi) > safe*err2
	  J(fi,xi) = deriv; 
	  mask(fi) = 0;
	end
      end
    end
    if (mask == 0) then
      break
    end
    j = j+1;
    if j == maxit
      error('DIFFEXT didn''t converge. Try to increase gstep_ (default 0.01)')
    end
    D1 = D2;
  end
  [m_err,i] = max(err);
  if m_err > 1e-12
    error('DIFFEXT obtains an accuracy > 1e-12. Try to increase gstep_ (default 0.01)')
  end
//  [m_err,i] = max(relerr);
//  if m_err > 1e-12
//    dyn_disp(err)
//    dyn_disp(relerr)
//    dyn_disp(D2) 
//    error('DIFFEXT obtains an relative accuracy > 1e-12. Try to increase gstep_ (default 0.01)')
//  end
end
 
// 10/12/2001 MJ modified initial h
// 02/25/2002 MJ put equation look inside
 
 
 
 
 
 
