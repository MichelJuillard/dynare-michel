var y,x;
varexo u,v;
parameters a, b, c, d, e, m, n;

a=-0.8;
b=0.9;
c=0.9;
d=1;
e=1;
m=50;
n=0.2;

model;
x=b*x(-1)+u;
a*y(+1)+y-(a*b^4+1)*x(-1)^4-4*a*b*x(-1)*u^3-4*a*b^3*x(-1)^3*u-6*a*(b*x(-1)*u)^2-a*u^4-v;
end;

initval;
x=0;
y=0;
u=0;
v=0;
end;

vcov=[1 0; 0 1];

order = 4;
