function [H]=hessext(f,x,options,varargin)
[nargout,nargin] = argn(0)
//---------------------------------------------------------------------------
//HESSEXT   Numerical approximation for hessian.
//          The method is Richardson`s extrapolation.
// Sample call
//   [H] = hessext(f,x,options,varargin)
// Inputs
//   f        name of the function
//   x        differentiation point
//   options  matrix of algorithm parameters
//   %delta    error goal (1e-12)
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
    %delta = options(1);
    toler = options(2);
    maxit = options(3)
  else
    %delta = 1e-9;
    toler = 1e-9;
    maxit = 20;
  end
else
  %delta = 1e-9;
  toler = 1e-9;
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
H = sparse([],[],[nf,nx*nx]);
for xi = 1:nx
  for xj = xi:nx
    err = big*ones(nf,1);
    relerr = big*ones(nf,1);
    err2 = big;
    hx = max(abs(x(xi))/5,gstep_);
    hy = max(abs(x(xj))/5,gstep_);
    dx = zeros(nx,1);
    dy = zeros(nx,1);
    dx(xi) = hx;
    dy(xj) = hy;
    if nargin > 3
      fss = evstr(f+'(x+dx+dy,varargin)');
      fsm = evstr(f+'(x+dx-dy,varargin)');
      fms = evstr(f+'(x-dx+dy,varargin)');
      fmm = evstr(f+'(x-dx-dy,varargin)');
    else	
      fss = evstr(f+'(x+dx+dy)');
      fsm = evstr(f+'(x+dx-dy)');
      fms = evstr(f+'(x-dx+dy)');
      fmm = evstr(f+'(x-dx-dy)');
    end
    D1 = (fss-fsm-fms+fmm)/(4*hx*hy);
    mask = ones(nf,1);
    j = 2;
    while (1)
      D2 = zeros(nf,j);
      hx = hx/con;
      hy = hy/con;
      dx = zeros(nx,1);
      dy = zeros(nx,1);
      dx(xi) = hx;
      dy(xj) = hy;
      if nargin > 3
	fss = evstr(f+'(x+dx+dy,varargin)');
	fsm = evstr(f+'(x+dx-dy,varargin)');
	fms = evstr(f+'(x-dx+dy,varargin)');
	fmm = evstr(f+'(x-dx-dy,varargin)');
      else	
	fss = evstr(f+'(x+dx+dy)');
	fsm = evstr(f+'(x+dx-dy)');
	fms = evstr(f+'(x-dx+dy)');
	fmm = evstr(f+'(x-dx-dy)');
      end
      for fi = 1:nf
	if mask(fi)
	  err2 = big;
	  D2(fi,1) = (fss(fi)-fsm(fi)-fms(fi)+fmm(fi))/(4*hx*hy);
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
	  if err(fi)  < toler | err(fi) > safe*err2
	    H(fi,(xi-1)*nx+xj) = deriv;
	    H(fi,(xj-1)*nx+xi) = deriv;
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
    if m_err > toler
      pause
      dyn_disp(D2)
      dyn_disp(err)
      dyn_disp([x dx dy])
      error('HESSEXT obtains an accuracy > 1e-12. Try to increase gstep_ (default 0.01)')
    end
  end
end

// 10/12/2001 MJ modified initial h
 
 
 
 
 
