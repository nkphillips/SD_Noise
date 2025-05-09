% estimateResponseBias  
% 
% This function estimates the response bias of the subject using the response bias model.
% 
% Inputs:
%   fixed_params: a cell array containing the fixed parameters
%   p: a struct containing the model parameters
% 
% Outputs:
%   response_bias: a struct containing the estimated parameters, SSE, and R2

function response_bias = estimateResponseBias(init_params, fixed_params, p)

    %% Grid search for starting parameters  
    
    start_params = gridSearchResponseBias(init_params, fixed_params, 'coarse', p.response_bias_bounds);
    [start_params, start_sse] = gridSearchResponseBias(start_params, fixed_params, 'fine', p.response_bias_bounds);
    free_params = start_params;

    %% Define the model function handle

    response_bias_model = @(free_params) calcResponseBiasFit(free_params, fixed_params);

    %% Fit the model

    [params_est, sse, exit_flag] = fmincon(response_bias_model, free_params, [], [], [], [], p.response_bias_bounds(2,:), p.response_bias_bounds(1,:), [], p.fmincon_options);

    %% Evaluate the model

    mu = params_est(1);
    sigma = params_est(2);

    est_p_CW = calc_pCW(fixed_params{1}(:,1), mu, sigma, fixed_params{2});
    r2 = calcR2(fixed_params{1}(:,2), est_p_CW);

    %% Compile response bias

    response_bias.start_params = start_params;
    response_bias.start_sse = start_sse;
    response_bias.params_est = params_est;
    response_bias.sse = sse;
    response_bias.r2 = r2;
    response_bias.exit_flag = exit_flag;
    
end 