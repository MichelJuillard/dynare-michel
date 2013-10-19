var y a k c i h eoy eoc eoh oy oc oh;
varexo e eeoy eeoc eeoh;

parameters theta rho eta gam bet delta aa r11 r12 r13 r21 r22 r23 r31 r32 r33 scy shc shy;

bet = 0.99;
delta = 0.025;
theta = 0.2;
rho = 0.9959;
eta = 1.0051;
gam = 0.0045;
aa = 1.8;
r11 = 0.99;
r12 = 0;
r13 = 0;
r21 = 0;
r22 = 0.99;
r23 = 0;
r31 = 0;
r32 = 0;
r33 = 0.99;
scy = 0.0040;
shy = 0.0015;
shc = 0.0010;

model(bytecode, block,cutoff=0);
exp(y) = exp(a)*exp(k(-1))^theta*exp(h)^(1-theta);
a = (1-rho)*aa+rho*a(-1)+e;
exp(y) = exp(c) + exp(i);
eta*exp(k) = (1-delta)*exp(k(-1))+exp(i);
gam*exp(c)*exp(h) = (1-theta)*exp(y);
eta/exp(c) = bet*(1/exp(c(+1)))*(theta*(exp(y(+1))/exp(k))+1-delta);
eoy = r11*eoy(-1) + r12*eoc(-1) + r13*eoh(-1) + eeoy;
eoc = r21*eoy(-1) + r22*eoc(-1) + r23*eoh(-1) + scy*eeoy+eeoc;
eoh = r31*eoy(-1) + r32*eoc(-1) + r33*eoh(-1) + shy*eeoy+shc*eeoc+eeoh;
oy = y + eoy;
oc = c + eoc;
oh = h + eoh;
end;

initval;
e=0;
eeoy=0;
eeoc=0;
eeoh=0;
y=    7.99331700544506;
a=    1.8;
k=    9.59646163090336;
c=    7.83132048725623;
i=    6.09323152367607;
h=    5.34253084908048;
eoy=   0.0000;
eoc=   0.0000;
eoh=         0;
oy=    7.99331700544506;
oc=    7.83132048725623;
oh=    5.34253084908048;
k=k+0.000001;

end;

options_.dynatol.f=1e-12;
options_.slowc=1;
steady(solve_algo=2,maxit=50);

options_.dynatol.f=4e-5;

shocks;
var e;
periods 1;
values 0.02;
end;

model_info;

simul(periods=2000, stack_solve_algo = 0,maxit=10);
rplot y;
rplot k;
