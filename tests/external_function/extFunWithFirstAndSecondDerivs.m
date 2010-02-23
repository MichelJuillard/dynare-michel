function [y df d2f]=extFunWithFirstAndSecondDerivs(a,b)
y=a*(b^2);

da=b^2;
db=2*a*b;
df=[da db];

d2f=[0 2*b; 2*b 2*a];
end