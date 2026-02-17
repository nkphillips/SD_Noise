% calcResponseBiasFit  
% 
% This function calculates the negative log-likelihood of the response bias model.
% 
% Inputs:
%   free_params: [mu, sigma]
%   fixed_params: [probe_offsets, responses]
%   guess_rate: the guess rate of the subject
% 
% Outputs:
%   nll: the negative log-likelihood of the response bias model

function nll = calcResponseBiasFit(free_params, fixed_params, guess_rate)

    %% Extract parameters 

    probe_offsets = fixed_params(:,1);
    responses = fixed_params(:,2);
    
    mu = free_params(1);
    sigma = free_params(2);

    %% Calculate pCW

    p_CW = calc_pCW(probe_offsets, mu, sigma, guess_rate);

    %% Calculate NLL

    nll = calcNLL(responses, p_CW);

end