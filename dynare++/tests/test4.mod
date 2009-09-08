var y,x;
varexo u,v;
parameters a, b, c, d, e, f, m;

a=0.8;
b=0.9;
c=0.9;
d=1;
e=1;
m=50;
f = 1;

model;
x = a*x(-1)+u;
c*y(+1)^2+d*y^2+e*x^2-(c+d)*m^2-(c*b*b*a*a+d*b*b+e*a*a)*x(-1)^2-(c*b*b+e)*u^2-2*(c*m*b*a+d*m*b)*x(-1)-2*c*m*b*u-2*(c*b*b*a+e*a)*x(-1)*u-d*f^2*v^2-2*d*m*f*v-2*d*b*f*x(-1)*v=0;
end;

initval;
x=1;
y=21;
u=0;
v=0;
end;

vcov=[1 0; 0 1];

order = 2;
