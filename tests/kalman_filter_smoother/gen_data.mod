var w x y z dw dx dy;
varexo e_w e_x e_y e_z;

parameters rho_w rho_x rho_y rho_z a1 a2 a3 b c;

rho_w = 0.5;
rho_x = 0.5;
rho_y = 0.5;
rho_z = 0.8;

a1 = 0.1;
a2 = -0.1;
a3 = 0.1;
b = 1;
c = 1;

model(linear);
dw = rho_w*dw(-1)+a1*(x(-1)-b*y(-1)-c*w(-1))+e_w;
dx = rho_x*dx(-1)+a2*(x(-1)-b*y(-1)-c*w(-1))+e_x;
dy = rho_y*dy(-1)+a3*(x(-1)-b*y(-1))+e_y;
z = rho_z*z(-1)+dw-dx+e_z;
w = w(-1)+dw;
x = x(-1)+dx;
y = y(-1)+dy;
end;

shocks;
var e_w; stderr 0.05;
var e_x; stderr 0.05;
var e_y; stderr 0.05;
var e_z; stderr 0.05;
end;

stoch_simul(periods=2000,irf=0);

plot([w x y z]);

save data.mat w x y z dw dx dy;