% computes the steady state of fs2000 analyticaly
% largely inspired by the program of F. Schorfheide
function [ys,check] = fs2000a_steadystate(junk,ys)
  global alp bet gam mst rho psi del;

  check = 0;
  
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  p  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/p;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  r  = mst/bet;
  w  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
  
  P_obs = 1;
  Y_obs = 1;
  
  ys =[
      c
      d
      dA
      e
      gp_obs
      gy_obs
      k
      l
      m
      n
      p
      P_obs
      r
      w
      y
      Y_obs
      ];