function [x, options] = fstop(funfcn,x,options,grad,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)
% FSTOP  Optimisation with user intervention

if nargin<3, options = []; end
options = foptions(options);
prnt = options(1);
tol = options(2);
tol2 = options(3);
evalstr = [funfcn];
stop_panel = P10;
if ~any(funfcn<48)
    evalstr=[evalstr, '(x'];
    for i=1:nargin - 4
        evalstr = [evalstr,',P',int2str(i)];
    end
    evalstr = [evalstr, ')'];
end
n = prod(size(x));
if (~options(14)) 
    options(14) = 200*n; 
end
xin = x(:);
v = xin;
x(:) = v; fv = eval(evalstr);
usual_delta = 0.05;
zero_term_delta = 0.00025;
for j = 1:n
   y = xin;
   if y(j) ~= 0
      y(j) = (1 + usual_delta)*y(j);
   else
      y(j) = zero_term_delta;
   end
   v = [v y];
   x(:) = y; f = eval(evalstr);
   fv = [fv  f];
end
[fv,j] = sort(fv);
v = v(:,j);
cnt = n+1;
if prnt
   clc
   format compact
   format short e
   home
   cnt
   disp('initial ')
   disp(' ')
   v
   f
end
alpha = 1;  beta = 1/2;  gamma = 2;
[n,np1] = size(v);
onesn = ones(1,n); 
ot = 2:n+1;
on = 1:n;
while cnt < options(14)
    if max(max(abs(v(:,ot)-v(:,onesn)))) <= tol & ...
           max(abs(fv(1)-fv(ot))) <= tol2
        break
    end
    vbar = (sum(v(:,on)')/n)';
    vr = (1 + alpha)*vbar - alpha*v(:,n+1);
    x(:) = vr;
    fr = eval(evalstr); 
    cnt = cnt + 1; 
    vk = vr;  fk = fr; how = 'reflect ';
    if fr < fv(n)
        if fr < fv(1)
            ve = gamma*vr + (1-gamma)*vbar;
            x(:) = ve;
            fe = eval(evalstr);
            cnt = cnt + 1;
            if fe < fv(1)
                vk = ve; fk = fe;
                how = 'expand  ';
            end
        end
    else
        vt = v(:,n+1); ft = fv(n+1);
        if fr < ft
            vt = vr; ft = fr;
        end
        vc = beta*vt + (1-beta)*vbar;
        x(:) = vc;
        fc = eval(evalstr); 
        cnt = cnt + 1;
        if fc < fv(n)
            vk = vc; fk = fc;
            how = 'contract';
        else
            for j = 2:n
                v(:,j) = (v(:,1) + v(:,j))/2;
                x(:) = v(:,j);
                fv(j) = eval(evalstr); 
            end
        cnt = cnt + n-1;
        vk = (v(:,1) + v(:,n+1))/2;
        x(:) = vk;
        fk = eval(evalstr); 
        cnt = cnt + 1;
        how = 'shrink  ';
        end
    end
    v(:,n+1) = vk;
    fv(n+1) = fk;
    [fv,j] = sort(fv);
    v = v(:,j);
    if prnt
        home
        cnt
        disp(how)
        disp(' ')
        v
        fv
    end
   if ~isempty(stop_panel)
     stop_panel=get(gcf,'UserData');
     if stop_panel==1;
       options(14)=cnt;
       disp(' ')
       disp(['Optimisation terminated by user after ', ...
                int2str(options(14)),' function evaluations.']);
     end
   end
end
x(:) = v(:,1);
if prnt, format, end
options(10)=cnt;
options(8)=min(fv); 
if cnt==options(14) 
    if options(1) >= 0
        if get(gcf,'UserData')~=1
          disp(['Warning: Maximum number of iterations (', ...
                 int2str(options(14)),') has been exceeded']);
          disp( '         (increase OPTIONS(14)).')
        end
    end
end

% end of m-file