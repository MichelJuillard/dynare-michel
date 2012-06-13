% computes the steady state of fs2000 analyticaly
% largely inspired by the program of F. Schorfheide
function [ys,check] = fs2000_b1L1L_steadystate(ys,exe)
  global M_
  
  alp = M_.params(1); 
  bet = M_.params(2); 
  gam = M_.params(3); 
  mst = M_.params(4); 
  rho = M_.params(5); 
  psi = M_.params(6); 
  del = M_.params(7); 

  check = 0;
  
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
  
  P_obs = 1;
  Y_obs = 1;
  
  P2=P;
  c2=c;
  
  ys =[
    m     
    P     
    c     
    e     
    W     
    R     
    k     
    d     
    n     
    l     
    gy_obs
    gp_obs
    Y_obs 
    P_obs 
    y     
    dA
    P2
    c2  ];