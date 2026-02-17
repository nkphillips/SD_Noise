function p_CW = calc_pCW(x, mu, sigma, guess_rate)
    % Handles both single and multiple parameter cases
    % 
    % Inputs:
    %   x: probe offsets (vector)
    %   mu: bias parameter(s) (scalar or vector)
    %   sigma: width parameter(s) (scalar or vector)
    %   guess_rate: guess rate (scalar)
    %
    % Outputs:
    %   p_CW: probability of CW response(s) (vector or matrix)

    if isscalar(mu) && isscalar(sigma)
        % Single parameter case (current usage)
        p_CW = (1 - guess_rate) * normcdf(x, mu, sigma) + (0.5 * guess_rate);
    else
        % Multiple parameter case (for vectorization)
        p_CW = zeros(length(x), length(mu));
        for i = 1:length(mu)
            p_CW(:,i) = (1 - guess_rate) * normcdf(x, mu(i), sigma(i)) + (0.5 * guess_rate);
        end
    end
    
    % Avoid log(0)
    p_CW = max(p_CW, 1e-10);
    
end 