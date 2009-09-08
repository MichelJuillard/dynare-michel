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
a*y(+1)+y-(a*b^3+1)*x(-1)^3-3*a*b*x(-1)*u^2-3*a*b^2*x(-1)^2*u-a*u^3-a-v^2;
end;

initval;
x=0;
y=0;
u=0;
v=0;
end;

vcov=[1 0; 0 1];

order = 3;
