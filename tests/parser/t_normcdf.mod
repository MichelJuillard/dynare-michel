var y1, y2, y3, x1, x2, x3;

model;
x1 = 1.96;
x2 = 1;
x3 = 0.5;
y1 = normcdf(x1(-1),0,1);
y2 = normcdf(-x1,-x2,1);
y3 = normcdf(x1/2,0,x3(+1));
end;

initval;
y1 = 0;
y2 = 0;
y3 = 0;
x1 = 0;
x2 = 0;
x3 = 1;
end;

steady;

if abs(oo_.steady_state(1) - pnorm(1.96,0,1)) > 1e-12;
   error('Error 1 in t_normcdf')
end;
if abs(oo_.steady_state(2) - pnorm(-1.96,-1,1)) > 1e-12;
   error('Error 2 in t_normcdf')
end;
if abs(oo_.steady_state(3) - pnorm(1.96/2,0,1/2)) > 1e-12;
   error('Error 3 in t_normcdf')
end;

z = [oo_.steady_state(4); oo_.steady_state; oo_.steady_state(6)];
[junk,JJ] = t_normcdf_dynamic(z,[]);

if abs(JJ(4,1) + dnorm(1.96,0,1)) > 1e-12;
   error('Error 4 in t_normcdf')
end;
if abs(JJ(5,5) - dnorm(-1.96,-1,1)) > 1e-12;
   error('Error 5 in t_normcdf')
end;
if abs(JJ(5,6) + dnorm(-1.96,-1,1)) > 1e-12;
   error('Error 6 in t_normcdf')
end;
if abs(JJ(6,5) + dnorm(1.96/2,0,1/2)/2) > 1e-12;
   error('Error 7 in t_normcdf')
end;
if abs(JJ(6,8) - (1/2)*((1.96/2)/(1/2)^2)*dnorm(1.96/2,0,1/2)) > 1e-12;
   error('Error 8 in t_normcdf')
end;

