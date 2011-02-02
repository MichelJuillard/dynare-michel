function [Ui,Vi,n0,np,ixmC0Pres] = exclusions(nvar,nexo,options_ms)

indxC0Pres = options_ms.cross_restrictions;
nlags = options_ms.nlags;

Qi = options_ms.Qi;
Ri1 = options_ms.Ri;

k = nlags*nvar+1;

Ri = zeros(k,k,nvar);
sR = size(Ri1);
Ri(1:sR(1),1:sR(2),1:sR(3)) = Ri1;

for n=1:nvar
 Ui{n} = null(Qi(:,:,n));
 Vi{n} = null(Ri(:,:,n));
 n0(n) = size(Ui{n},2);
 np(n) = size(Vi{n},2);
end

ixmC0Pres = NaN;