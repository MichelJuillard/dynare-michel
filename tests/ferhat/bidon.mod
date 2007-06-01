var c i g infl y k l r p1 q1 p2 q2 r0;
varexo g_bar;

parameters a b d e f h j m n o p q;
a=0.4*0.6;
b=0.3*0.6;
d=0.1;
e=0.15;
f=1;
h=0.15;
j=1;
m=1;
n=1;
o=1;

model(SPARSE_DLL,gcc_compiler,cutoff=1e-12);
/*0*/ k=(1-h)*k(-1)+i;  /*k:0*/
/*1*/ y=l^j*k^m;          /*l:1*/
/*2*/ c=y*a+b+0.3*c(-1)+0.1*c(+1)+0.*g_bar(-10);          /*c:2*/
/*3*/ infl=0.02*y+0.5*r;         /*infl:3*/
/*4*/ i=d*(y-y(-1))+e/**r*/;  /*i4*/
/*5*/ g=f*g_bar;              /*g:5*/
/*6*/ y=0.6*(c+i+g)+/*0.1*y(-2)+0.1*y(+2)+*/0.1*y(-1)+0.1*y(+1);          /*y:6*/
/*7*/ r=y-1+infl-0.02;         /*r7*/
/*8*/ p1=i+0.5*q1;
/*9*/ q1=0.5*p1+c;
/*10*/ q2=0.5*p2+r0;
/*11*/ p2=0.5*q2+p1;
/*12*/ r0=r;
end;

initval;
g_bar=0.15;
c=0.7;
i=0.15;
g=0.15;
y=1;
k=1;
l=1;
infl=0.02;
r=0;
r0=r;
p1=2/3;
q1=3.1/3;
q2=4/9;
p2=8/9;

end;

steady(solve_algo=2);
//check;

shocks;
var g_bar;
periods 1;
values 0.16;
end;

options_.slowc = 1;


simul(periods=80);


rplot c;
rplot y;
