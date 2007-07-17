var dx dy;
varobs dx dy;

bvar_density(datafile = test, first_obs = 20, bvar_prior_flat,
             bvar_prior_train = 10) 8;

bvar_forecast(forecast = 10, bvar_replic = 10000, nobs = 200) 8;
