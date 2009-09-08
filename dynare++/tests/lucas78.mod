var Y P;

varexo EXO_Y;

parameters beta gamma rho y_ss;

beta = 0.95;
gamma= 0.5;
rho  = 0.9;
y_ss = 2;

model;
Y-y_ss = rho*(Y(-1)-y_ss) + EXO_Y;
Y^(-gamma)*P = beta*Y(+1)^(-gamma)*(P(+1) + Y(+1));
end;

initval;
Y = 2;
P = 38;
end;

vcov = [
10
];

order = 7;
