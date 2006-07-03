function [x,f,fvec,%check]=lnsrch(xold,fold,g,p,stpmax,func,varargin)
// Copyright (C) 2001 Michel Juillard
// 
 
alf = .0001;
tolx = .000000000037;
alam = 1;
 
nn = size(xold,1);
 
summ = sqrt(sum(p .* p));
if abs(summ) > abs(stpmax) then
  p = p .* stpmax/summ;
end
 
slope = g'*p;
 
test = max(abs(p)' ./ max([abs(xold)';ones(1,nn)],'r'));
alamin = tolx/test;
 
if alamin>.1 then
  alamin = .1;
end
 
while 1 then
  if alam<alamin then
    %check = 1;
    break
     ;
     
  end
   
  x = xold+alam*p;

  if argn(2) > 6
    fvec = evstr(func+'(x,varargin)');
  else
    fvec = evstr(func+'(x)');
  end

  f = real(.5*fvec'*fvec);
  if or(isnan(fvec)) then
    alam = alam/2;
    alam2 = alam;
    f2 = f;
    fold2 = fold;
  else
    if abs(f) <= abs(fold+alf*alam*slope) then
      %check = 0;
      break;
    else
      if alam==1 then
        tmplam = -slope/(2*(f-fold-slope));
      else
        rhs1 = f-fold-alam*slope;
        rhs2 = f2-fold2-alam2*slope;
        a = (rhs1/(alam^2)-rhs2/(alam2^2))/(alam-alam2);
        b = (-alam2*rhs1/(alam^2)+alam*rhs2/(alam2^2))/(alam-alam2);
        if a==0 then
          tmplam = -slope/(2*b);
        else
          disc = b^2-3*a*slope;
           
          if abs(disc)<0 then
            error('Roundoff problem in nlsearch');
          else
            tmplam = (-b+sqrt(disc))/(3*a);
          end
           
        end
         
        if abs(tmplam) >.5*abs(alam) then
          tmplam = .5*alam;
        end
         
      end
       
      alam2 = alam;
      f2 = f;
      fold2 = fold;
      alam = max([abs(tmplam);.1*abs(alam)]);
       
    end
  end
end
 
// 01/14/01 MJ lnsearch is now a separate function
