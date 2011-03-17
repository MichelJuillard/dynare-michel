function [yy, xdir, isig, lam]=log_trans_(y0,xdir0)

if nargin==1,
  xdir0='';
end
f=inline('skewness(log(y+lam))','lam','y');
isig=1;
if ~(max(y0)<0 | min(y0)>0)
  if skewness(y0)<0,
    isig=-1;
    y0=-y0;
  end
  n=hist(y0,10);
  if n(1)>20*n(end),
    try lam=fzero(f,[-min(y0)+10*eps -min(y0)+abs(median(y0))],[],y0);
    catch
      yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
      if abs(yl(1))<abs(yl(2))
        lam=-min(y0)+eps;
      else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
      end
    end
    yy = log(y0+lam);
    xdir=[xdir0,'_logskew'];
  else
    isig=0;
    lam=0;
    yy = log(y0.^2);
    xdir=[xdir0,'_logsquared'];
  end
else
  if max(y0)<0
    isig=-1;    
    y0=-y0;
    %yy=log(-y0);
    xdir=[xdir0,'_minuslog'];
  elseif min(y0)>0
    %yy=log(y0);
    xdir=[xdir0,'_log'];
  end
  try lam=fzero(f,[-min(y0)+10*eps -min(y0)+median(y0)],[],y0);
  catch
    yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
    if abs(yl(1))<abs(yl(2))
      lam=-min(y0)+eps;
    else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
    end
  end
  yy = log(y0+lam);
end
