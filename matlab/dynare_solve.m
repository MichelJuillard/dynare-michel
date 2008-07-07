function [x,info] = dynare_solve(func,x,jacobian_flag,varargin)

% function [x,info] = dynare_solve(func,x,jacobian_flag,varargin)
% proposes different solvers
%
% INPUTS
%    func:             name of the function to be solved
%    x:                guess values
%    jacobian_flag=1:  jacobian given by the 'func' function
%    jacobian_flag=0:  jacobian obtained numerically
%    varargin:         list of arguments following jacobian_flag
%    
% OUTPUTS
%    x:                solution
%    info=1:           the model can not be solved
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2001-2008)
% Gnu Public License.


  global options_
  
  options_ = set_default_option(options_,'solve_algo',2);
  info = 0;
  if options_.solve_algo == 0
    if ~isempty(which('fsolve'))
      options=optimset('fsolve');
      options.MaxFunEvals = 50000;
      options.MaxIter = 2000;
      options.TolFun=1e-8;
      options.Display = 'iter';
      if jacobian_flag
	options.Jacobian = 'on';
      else
	options.Jacobian = 'off';
      end
      [x,fval,exitval,output] = fsolve(func,x,options,varargin{:});
      if exitval > 0
	info = 0;
      else
	info = 1;
      end
      return
    else 
      options_.solve_algo = 1;
    end
  end

  if options_.solve_algo == 1
    nn = size(x,1);
    [x,info]=solve1(func,x,1:nn,1:nn,jacobian_flag,varargin{:});
  elseif options_.solve_algo == 2
    nn = size(x,1) ;
    tolf = options_.solve_tolf ;

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
    
%    f = 0.5*fvec'*fvec ;

    if max(abs(fvec)) < tolf
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
      [x,info]=solve1(func,x,j1(r(i):r(i+1)-1),j2(r(i):r(i+1)-1),jacobian_flag,varargin{:});
      if info & options_.debug
	error(sprintf('Solve block = %d check = %d\n',i,info));
      end
    end
    fvec = feval(func,x,varargin{:});
    if max(abs(fvec)) > tolf
      [x,info]=solve1(func,x,1:nn,1:nn,jacobian_flag,varargin{:});
    end
  elseif options_.solve_algo == 3
    if jacobian_flag
      [x,info] = csolve(func,x,func,1e-6,500,varargin{:});
    else
      [x,info] = csolve(func,x,[],1e-6,500,varargin{:});
    end 
  end
%    fvec1 = feval(func,x,varargin{:})

  % 08/28/03 MJ add a final call to solve1 for solve_algo == 1 in case
  %             initvals generates 'false' zeros in the Jacobian
  