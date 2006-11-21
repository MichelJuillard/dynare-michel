function [xcum, paracdf] = cumdens(para, pshape, p1, p2, p3, p4)
% This procedure transforms x vectors into cumulative values 
% pshape: 0 is point mass, both para and p2 are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu)
%         5 is UNIFORM [p1,p2]
% Adapted by M. Ratto from MJ priordens.m

nprio 	= length(pshape);

i = 1;
while i <=  nprio;
	a = 0;	
	b = 0;
	if pshape(i) == 1;     % (generalized) BETA Prior 
		mu = (p1(i)-p3(i))/(p4(i)-p3(i));
		stdd = p2(i)/(p4(i)-p3(i));
		a = (1-mu)*mu^2/stdd^2 - mu;
		b = a*(1/mu - 1);
		%lnprior = lnprior + lpdfgbeta(para(i),a,b,p3(i),p4(i))   ;
		para(:,i) = (para(:,i)-p3(i))./(p4(i)-p3(i));
		xcum(:,i) = betacdf(para(:,i),a,b)   ;
    as=num2str(a,'%20.16e');
    bs=num2str(b,'%20.16e');
    p3s=num2str(p3(i),'%20.16e');
    p4s=num2str(p4(i),'%20.16e');
    paracdf{i} = inline(['betacdf((X-',p3s,')./(',p4s,'-',p3s,'),',as,',',bs,')'],'X');
	elseif pshape(i) == 2; % GAMMA PRIOR 
     	b = p2(i)^2/(p1(i)-p3(i));
		a = (p1(i)-p3(i))/b;
		%lnprior = lnprior + lpdfgam(para(i)-p3(i),a,b);
		xcum(:,i) = gamcdf(para(:,i)-p3(i),a,b);
    as=num2str(a,'%20.16e');
    bs=num2str(b,'%20.16e');
    p3s=num2str(p3(i),'%20.16e');
    paracdf{i} = inline(['gamcdf((X-',p3s,'),',as,',',bs,')'],'X');
	elseif pshape(i) == 3; % GAUSSIAN PRIOR 
     %lnprior = lnprior + lpdfnorm(para(i),p1(i),p2(i));
     xcum(:,i) = normcdf(para(:,i),p1(i),p2(i));
    p1s=num2str(p1(i),'%20.16e');
    p2s=num2str(p2(i),'%20.16e');
    paracdf{i} = inline(['normcdf(X,',p1s,',',p2s,')'],'X');
	elseif pshape(i) == 4; % INVGAMMA1 PRIOR 
     	%lnprior = lnprior + lpdfig1(para(i),p1(i),p2(i));
      xcum(:,i)=para(:,i);
	elseif pshape(i) == 5; % UNIFORM PRIOR 
     	%lnprior = lnprior + log(1/(p2(i)-p1(i)));
  		xcum(:,i) = (para(:,i)-p1(i))./(p2(i)-p1(i));
    p1s=num2str(p1(i),'%20.16e');
    p2s=num2str(p2(i),'%20.16e');
    paracdf{i} = inline(['(X-',p1s,')./(',p2s,'-',p1s,')'],'X');
% 	elseif pshape(i) == 6; % INVGAMMA2 PRIOR 
%     	lnprior = lnprior + lpdfig2(para(i),p1(i),p2(i));
	end;
	i = i+1;
end;

