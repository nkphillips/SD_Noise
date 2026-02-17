% calcSerialDependenceFit
% 
% Calculates metric (SSE or NLL) for serial dependence model
% 
% Inputs:
%   free_params: [amplitude, width, baseline]
%   fixed_params: [probe_offsets, responses, delta_thetas]
%   p: parameter struct (expects p.sd_objective and p.guess_rate)
% 
% Outputs:
%   metric: the sum of squared errors or negative log-likelihood

function metric = calcSerialDependenceFit(free_params, fixed_params, p)

    %% Extract parameters
    
    amplitude = free_params(1);
    width = free_params(2);
    baseline = free_params(3);
    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'nll')
        sigma = free_params(4);
    else
        sigma = [];
    end

    %% Extract trial data
    
    probe_offsets = fixed_params(:,1);
    responses = fixed_params(:,2);
    delta_thetas = fixed_params(:,3);

    %% Calculate predicted bias using DoG
    
    predicted_bias = calcDoG(delta_thetas, [amplitude, width, baseline]);

    %% Calculate metric

    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
        metric = calcSSE(responses, predicted_bias);
    else
        p_CW = calc_pCW(probe_offsets, predicted_bias, sigma, p.guess_rate);
        metric = calcNLL(responses, p_CW);
    end
    
end 