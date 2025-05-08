% estimateResponseBias
    

function response_bias = estimateResponseBias(fixed_params, p, toggles)

    probe_offsets = fixed_params{1}(:,1);
    p_CW = fixed_params{1}(:,2);
    guess_rate = fixed_params{2};

    %% Grid search for starting parameters  
    
    start_params = grid_search(probe_offsets, p_CW, guess_rate);
    free_params = start_params;

    %% Define the model function handle

    response_model = @(free_params) calc_fit(free_params, fixed_params);

    %% Fit the model

    [params_est, sse, exit_flag] = fmincon(response_model, free_params, [], [], [], [], p.response_bias_bounds(2,:), p.response_bias_bounds(1,:), [], p.fmincon_options);

    %% Evaluate the model

    mu = params_est(1);
    sigma = params_est(2);

    est_pCW = calc_pCW(probe_offsets, mu, sigma, guess_rate);
    r2  = 1 - sum((p_CW - est_pCW).^2) / sum((p_CW - mean(p_CW)).^2);

    %% Compile response bias

    response_bias.start_params = start_params;
    response_bias.start_sse = sse;
    response_bias.start_r2 = r2;
    response_bias.params_est = params_est;
    response_bias.sse = sse;
    response_bias.r2 = r2;
    response_bias.exit_flag = exit_flag;
    
end 