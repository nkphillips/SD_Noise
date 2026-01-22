% estimateResponseBias  
% 
% This function estimates the response bias of the subject using the response bias model.
% 
% Inputs:
%   fixed_params: [probe_offsets, responses]
%   p: a struct containing the model options, parameter bounds
% 
% Outputs:
%   start_params, start_nll, params_est, nll, exit_flag

function [start_params, start_nll, params_est, nll, exit_flag] = estimateResponseBias(fixed_params, p)

    %% Grid search for starting parameters

    start_params = gridSearchResponseBias(p.rb_init_params, fixed_params, 'coarse', p);
    [start_params, start_nll] = gridSearchResponseBias(start_params, fixed_params, 'fine', p);
    free_params = start_params;

    %% Define the model function handle

    rb_model = @(free_params) calcResponseBiasFit(free_params, fixed_params, p.guess_rate);

    %% Fit the model

    [params_est, nll, exit_flag] = fmincon(rb_model, free_params, [], [], [], [], p.rb_bounds(2,:), p.rb_bounds(1,:), [], p.fmincon_options);

end 