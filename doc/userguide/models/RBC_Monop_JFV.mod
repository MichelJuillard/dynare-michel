
// Adapted from Jesus Fernandez-Villaverde, Basic RBC Model with Monopolistic Competion Philadelphia, March 3, 2005






var y c k i l y_l w r  z;


varexo e;





parameters beta psi delta alpha rho gamma sigma epsilon;



alpha   = 0.33;


beta    = 0.99;


delta   = 0.023;


psi     = 1.75;


rho     = 0.95;  


sigma   = (0.007/(1-alpha));


epsilon = 10;





model;


  (1/c) = beta*(1/c(+1))*(1+r(+1)-delta);


  psi*c/(1-l) = w;


  c+i = y;


  y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha);


  w = y*((epsilon-1)/epsilon)*(1-alpha)/l;


  r = y*((epsilon-1)/epsilon)*alpha/k(-1);


  i = k-(1-delta)*k(-1);


  y_l = y/l;


  z = rho*z(-1)+e;


end;



initval;


  k = 9;


  c = 0.76;


  l = 0.3;


  w = 2.07;


  r = 0.03;


  z = 0; 


  e = 0;


end;

steady;

check;




shocks;


var e = sigma^2;


end;





steady;





stoch_simul(periods=2100);


