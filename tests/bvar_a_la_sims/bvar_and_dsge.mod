var dx dy;
varexo e_x e_y;

parameters rho_x rho_y;

rho_x = 0.5;
rho_y = -0.3;

model;
dx = rho_x*dx(-1)+e_x;
dy = rho_y*dy(-1)+e_y;
end;

estimated_params;
rho_x,NORMAL_PDF,0.5,0.1;
rho_y,NORMAL_PDF,-0.3,0.1;
stderr e_x,INV_GAMMA_PDF,0.01,inf;
stderr e_y,INV_GAMMA_PDF,0.01,inf;
end;

varobs dx dy;

check;

estimation(datafile = bvar_sample, mh_replic = 1200, mh_jscale = 1.3,
           first_obs = 20);

bvar_density(bvar_prior_train = 10) 8;

bvar_forecast(forecast = 10, bvar_replic = 2000, nobs = 200) 8;
