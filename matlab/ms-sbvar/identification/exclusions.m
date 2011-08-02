function [Ui,Vi,n0,np,ixmC0Pres] = exclusions(nvar,nexo,options_ms)

indxC0Pres = options_ms.cross_restrictions;
nlags = options_ms.nlags;

Qi1 = options_ms.Qi;
Ri1 = options_ms.Ri;

Ui = cell(nvar,1);
Vi = cell(nvar,1);
n0 = zeros(nvar,1);
np = zeros(nvar,1);

k = nlags*nvar+1;

for n=1:nvar
 Qi{n} = zeros(nvar,nvar);
 sQ = size(Qi1{n});
 if all(sQ) > 0
     Qi{n}(1:sQ(1),1:sQ(2)) = Qi1{n};
 end
 Ri{n} = zeros(k,k);
 sR = size(Ri1{n});
 if all(sR) > 0
     Ri{n}(1:sR(1),1:sR(2)) = Ri1{n};
 end

 if options_ms.constants_exclusion
        Ri{n}(sR(1)+1,k) = 1;
 end

 Ui{n} = null(Qi{n});
 Vi{n} = null(Ri{n});
 n0(n) = size(Ui{n},2);
 np(n) = size(Vi{n},2);
end

ixmC0Pres = NaN;