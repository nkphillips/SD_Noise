% estimateSerialDependence
%
% Serial dependence estimation with smart initialization
%
% Inputs:
%   raw_data: [probe_offsets, responses, delta_thetas]
%   p: parameter structure
%   init_params: initial parameters
%
% Outputs:
%   start_params: starting parameters from grid search
%   start_metric: negative log-likelihood or SSE at start
%   params_est: estimated parameters
%   final_metric: final negative log-likelihood or SSE
%   exit_flag: optimization exit flag

function [start_params, start_metric, params_est, final_metric, exit_flag] = estimateSerialDependence(raw_data, p, init_params)

    %% Grid search for starting parameters

    start_params = gridSearchSerialDependence(init_params, raw_data, 'coarse', p);
    [start_params, start_metric] = gridSearchSerialDependence(start_params, raw_data, 'fine', p);
    free_params = start_params;

    %% Unified model function handle based on objective (single entry point)
    
    sd_model = @(free_params) calcSerialDependenceFit(free_params, raw_data, p);
    
    %% Fit the model

    % Ensure bounds match the parameter count implied by objective
    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
        % SSE: [A, w, b]
        ub = p.sd_bounds(1, 1:3);
        lb = p.sd_bounds(2, 1:3);
    else
        % NLL: [A, w, b, sigma]
        ub = p.sd_bounds(1, :);
        lb = p.sd_bounds(2, :);
    end
    
    [params_est, final_metric, exit_flag] = fmincon(sd_model, free_params, [], [], [], [], lb, ub, [], p.fmincon_options);

end 