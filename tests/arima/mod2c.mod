var dx dy coint_err;
varexo e_x e_y;

parameters rho_x rho_y b a1 a2;

rho_x = 0.5;
rho_y = -0.3;
b = 1;
a1 = -0.1;
a2 = 0.2;

model;
dx = rho_x*dx(-1)+a1*coint_err(-1)+e_x;
dy = rho_y*dy(-1)+a2*coint_err(-1)+e_y;
coint_err = dx-b*dy+coint_err(-1);
end;

estimated_params;
rho_x,NORMAL_PDF,0.5,0.1;
rho_y,NORMAL_PDF,-0.3,0.1;
b,NORMAL_PDF,1,0.1;
a1,NORMAL_PDF,-0.1,0.1;
a2,NORMAL_PDF,0.2,0.1;

stderr e_x,INV_GAMMA_PDF,0.01,inf;
stderr e_y,INV_GAMMA_PDF,0.01,inf;
end;

varobs dx dy;

estimation(datafile=data2,nobs=100,mh_replic=0);
