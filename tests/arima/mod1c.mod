var dx dy x y;
varexo e_x e_y;

parameters rho_x rho_y;

rho_x = 0.5;
rho_y = -0.3;

model;
dx = rho_x*dx(-1)+e_x;
dy = rho_y*dy(-1)+e_y;
x = x(-1)+dx;
y = y(-1)+dy;
end;

estimated_params;
rho_x,NORMAL_PDF,0.5,0.1;
rho_y,NORMAL_PDF,-0.3,0.1;
stderr e_x,INV_GAMMA_PDF,0.01,inf;
stderr e_y,INV_GAMMA_PDF,0.01,inf;
stderr x,INV_GAMMA_PDF,0.01,inf;
stderr y,INV_GAMMA_PDF,0.01,inf;
end;

varobs x y;
options_.unit_root_vars = {'x'; 'y'};
estimation(datafile=data1,nobs=1000,mh_replic=2000,lik_init=2);