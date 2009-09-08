var y,x;
varexo u,v;
parameters a, b, c, d, e, f, g, h, j;

a=0.8;
b=0.9;
c=0.9;
d=1;
e=-0.556875;
f=-0.172125;
g=-0.9;
h=-0.2754;
j=-1.8;


model;
x=a*x(-1)+u;
c*y(+1)^2+d*y^2+e*x^2+f*u^2-d*v^2+g+h*x(-1)*u+j*x(-1)*v=0;
end;

initval;
x=0;
y=0.7237469;
u=0;
v=0;
end;

vcov=[1 0; 0 1];

order = 2;