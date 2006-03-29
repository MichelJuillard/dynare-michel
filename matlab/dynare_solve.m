% Copyright (C) 2001 Michel Juillard
%
function [x,cheik] = dynare_solve(func,x,jacobian_flag,varargin)
  global options_
  
  options_ = set_default_option(options_,'solve_algo',2);
  cheik = 0;
  if options_.solve_algo == 0
    if ~isempty(which('fsolve')) & sscanf(version('-release'),'%d') >= 13;
      options=optimset('fsolve');
      options.MaxFunEvals = 20000;
      options.TolFun=1e-8;
      options.Display = 'off';
      if jacobian_flag
	options.Jacobian = 'on';
      else
	options.Jacobian = 'on';
      end
      [x,fval,exitval,output] = fsolve(func,x,options,varargin{:});
      if exitval > 0
	cheik = 0;
      else
	cheik = 1;
      end
      return
    else 
      options_.solve_algo = 1;
    end
  end

  if options_.solve_algo == 1
    nn = size(x,1) ;
    [x,cheik]=solve1(func,x,1:nn,1:nn,jacobian_flag,varargin{:});
  elseif options_.solve_algo == 2
    nn = size(x,1) ;
    tolf = eps^(2/3) ;

    if jacobian_flag
      [fvec,fjac] = feval(func,x,varargin{:});
    else
      fvec = feval(func,x,varargin{:});
      fjac = zeros(nn,nn) ;
    end

    i = find(~isfinite(fvec));
    
    if ~isempty(i)
      disp(['STEADY:  numerical initial values incompatible with the following' ...
	    ' equations'])
      disp(i')
      error('exiting ...')
    end
    
    f = 0.5*fvec'*fvec ;

    if max(abs(fvec)) < 0.01*tolf
      return ;
    end

    if ~jacobian_flag
      fjac = zeros(nn,nn) ;
      dh = max(abs(x),options_.gstep*ones(nn,1))*eps^(1/3);
      for j = 1:nn
	xdh = x ;
	xdh(j) = xdh(j)+dh(j) ;
	fjac(:,j) = (feval(func,xdh,varargin{:}) - fvec)./dh(j) ;
      end
    end

    [j1,j2,r,s] = dmperm(fjac);
    
    for i=length(r)-1:-1:1
      [x,cheik]=solve1(func,x,j1(r(i):r(i+1)-1),j2(r(i):r(i+1)-1),jacobian_flag,varargin{:});
      if cheik
	error(sprintf('Solve block = %d check = %d\n',i,cheik));
      end
    end
    [x,cheik]=solve1(func,x,1:nn,1:nn,jacobian_flag,varargin{:});
      
  end
%    fvec1 = feval(func,x,varargin{:})

  % 08/28/03 MJ add a final call to solve1 for solve_algo == 1 in case
  %             initvals generates 'false' zeros in the Jacobian
  