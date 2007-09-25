function bvar_forecast(nlags)

    global options_
    
    options_ = set_default_option(options_, 'bvar_replic', 2000);
    if options_.forecast == 0
        error('bvar_forecast: you must specify "forecast" option')
    end
    [ny, nx, posterior, prior, forecast_data] = bvar_toolbox(nlags);

        
    sims_no_shock = NaN(options_.forecast, ny, options_.bvar_replic);
    sims_with_shocks = NaN(options_.forecast, ny, options_.bvar_replic);
    
    S_inv_upper_chol = chol(inv(posterior.S));

    % Option 'lower' of chol() not available in old versions of
    % Matlab, so using transpose
    XXi_lower_chol = chol(posterior.XXi)';
    
    k = ny*nlags+nx;
    
    % Declaration of the companion matrix:
    Companion_matrix = diag(ones(ny*(nlags-1),1),-ny);
    
    p = 0;
    d = 0;
    while d<=options_.bvar_replic
        
        Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);

        % Option 'lower' of chol() not available in old versions of
        % Matlab, so using transpose
        Sigma_lower_chol = chol(Sigma)';

        Phi = rand_matrix_normal(k, ny, posterior.PhiHat, XXi_lower_chol, Sigma_lower_chol);
        
        % All the eigenvalues of the companion matrix have to be on or inside the unit circle
        Companion_matrix(1:ny,:) = Phi(1:ny*nlags,:)'; 
        test = (abs(eig(Companion_matrix)));
        if any(test>1.0000000000001)
            p = p+1;
            continue
        end
        d = d+1;
        
        % Without shocks
        lags_data = forecast_data.initval;
        
        for t = 1:options_.forecast
            X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.xdata(t, :) ];
            y = X * Phi;
            lags_data(1:end-1,:) = lags_data(2:end, :);
            lags_data(end,:) = y;
            sims_no_shock(t, :, d) = y;
        end
    
        % With shocks
        lags_data = forecast_data.initval;
        for t = 1:options_.forecast
            X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.xdata(t, :) ];
            shock = (Sigma_lower_chol * randn(ny, 1))';
            y = X * Phi + shock;
            lags_data(1:end-1,:) = lags_data(2:end, :);
            lags_data(end,:) = y;
            sims_with_shocks(t, :, d) = y;
        end
    end
    
    disp('')
    disp(['Some of the VAR models sampled from the posterior distribution'])
    disp(['were found to be explosive (' int2str(p) ').'])
    disp('')
    
    % Plot graphs
    sims_no_shock_mean = mean(sims_no_shock, 3);

    sort_idx = round((0.5 + [-options_.conf_sig, options_.conf_sig]/2) * options_.bvar_replic);
    
    sims_no_shock_sort = sort(sims_no_shock, 3);
    sims_no_shock_down_conf = sims_no_shock_sort(:, :, sort_idx(1));
    sims_no_shock_up_conf = sims_no_shock_sort(:, :, sort_idx(2));

    sims_with_shocks_sort = sort(sims_with_shocks, 3);
    sims_with_shocks_down_conf = sims_with_shocks_sort(:, :, sort_idx(1));
    sims_with_shocks_up_conf = sims_with_shocks_sort(:, :, sort_idx(2));

    dynare_graph_init(sprintf('BVAR forecasts (nlags = %d)', nlags), ny, {'b-' 'g-' 'g-' 'r-' 'r-'});
    
    for i = 1:ny
        dynare_graph([ sims_no_shock_mean(:, i) ...
                       sims_no_shock_up_conf(:, i) sims_no_shock_down_conf(:, i) ...
                       sims_with_shocks_up_conf(:, i) sims_with_shocks_down_conf(:, i) ], ...
                     options_.varobs(i, :));
    end
    
    dynare_graph_close;

    
    % Compute RMSE
    
    if ~isempty(forecast_data.realized_val)
        
        sq_err_cumul = zeros(1, ny);
        
        lags_data = forecast_data.initval;
        for t = 1:size(forecast_data.realized_val, 1)
            X = [ reshape(flipdim(lags_data, 1)', 1, ny*nlags) forecast_data.realized_xdata(t, :) ];
            y = X * posterior.PhiHat;
            lags_data(1:end-1,:) = lags_data(2:end, :);
            lags_data(end,:) = y;
            sq_err_cumul = sq_err_cumul + (y - forecast_data.realized_val(t, :)) .^ 2;
        end
        
        rmse = sqrt(sq_err_cumul / size(sq_err_cumul, 1));
        
        fprintf('RMSE of BVAR(%d):\n', nlags);
        
        for i = 1:size(options_.varobs, 1)
            fprintf('%s: %10.4f\n', options_.varobs(i, :), rmse(i));
        end
    end