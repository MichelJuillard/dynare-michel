function y = rndprior(bayestopt_)

pshape=bayestopt_.pshape;
pmean=bayestopt_.pmean;
p1=bayestopt_.p1;
p2=bayestopt_.p2;
p3=bayestopt_.p3;
p4=bayestopt_.p4;
 
for i=1:length(pmean),
    
    switch pshape(i)
        
     case 1 %'beta'
      mu = (pmean(i)-p3(i))/(p4(i)-p3(i));
      stdd = p2(i)/(p4(i)-p3(i));
      A = (1-mu)*mu^2/stdd^2 - mu;
      B = A*(1/mu - 1);
      y(1,i) = beta_rnd(1, A, B);
      y(1,i) = y(1,i) * (p4(i)-p3(i)) + p3(i);
      
     case 2 %'gamma'
      mu = pmean(i)-p3(i);
      B = mu/p2(i)^2;              %gamm_rnd uses 1/B instead of B as param.
      A = mu*B;
      y(1,i) = gamm_rnd(1, A, B) + p3(i);
      
     case 3 %'normal'
      MU = pmean(i);
      SIGMA = p2(i);
      y(1,i) = randn*SIGMA+ MU;
      
     case 4 %'invgamma'
      nu = p2(i);
      s = p1(i);
      y(1,i) = 1/sqrt(gamm_rnd(1, nu/2, s/2));    %gamm_rnd uses 1/B
                                                  %instead of B as param.
      
     case 5 %'uniform'
      y(1,i) = rand*(p2(i)-p1(i)) + p1(i);
      
    end
end

% initial version by Marco Ratto

