function [ys,check] = wlash1_old_ss_steadystate(ys0,exo)
    global M_
    
    check = 0;
    
    params = M_.params;
    alpha = params(1);
    beta = params(2);
    delta = params(3);
    gamm = params(4);
    phi1 = params(5);
    eta = params(6);
    a = params(7);
    b = params(8);
    rho = params(9);
    phi2 = params(10);
    Psi = params(11);
    thetass = params(12);
 
    pi = thetass-1;
    en = 1/3;
    eR = 1/beta;
    y_k = (1/alpha)*(1/beta-1+delta);
    ek = en*y_k^(-1/(1-alpha));
    ec = ek*(y_k-delta);
    em = ec*(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^(-1/b);
    ey = ek*y_k;
    Xss = a*ec^(1-b)*(1+(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^((b-1)/b));
    Psi = (1-alpha)*(ey/en)*Xss^((b-phi1)/(1-b))*a*ec^(-b)*(1-en)^eta;
    n = log(en);
    k = log(ek);
    m = log(em);
    c = log(ec);
    y = log(ey);
    R = log(eR);
    z = 0;
    u = 0;
    
    ys = [y c k m n R pi z u]';
    M_.params(11) = Psi;
    
