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

model;
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
a = 1.7;
y = 8;
c = 8;
k = 10;
i = 5;
h = 4;
eoy = 0;
eoc = 0;
eoh = 0;
oy = y;
oc = c;
oh = h;
end;

//steady;
//check;

estimated_params;
theta , 0.22, 0.1, 0.5;
rho , 0.99, 0.7, 0.9999; 
eta , 1.0051, 1, 1.03;  
gam , 0.0045, 0.001, 0.01;
aa , 1.8, 0.1, 4;
r11 , 1.4187, -2, 2;
r12 , 0.2251, -2, 2;
r13 , -0.4441, -2, 2;
r21 , 0.0935, -2, 2;
r22 , 1.0236, -2, 2;
r23 , -0.0908, -2, 2;
r31 , 0.7775, -2, 2;
r32 , 0.3706, -2, 2;
r33 , 0.2398, -2, 2;
scy , 0.0040, -2, 2;
shy , 0.0015, -2, 2;
shc , 0.0010, -2, 2;
stderr e , 0.0056, 0, 0.2;
stderr eeoy , 0.0070, 0, 0.1;
stderr eeoc , 0.0069, 0, 0.1;
stderr eeoh , 0.0018, 0, 0.1;
end;

varobs oy oc oh;

observation_trends;
oy (log(eta));
oc (log(eta));
end;

//options_.debug=1;
estimation(datafile=idata,mode_compute=4,nograph);