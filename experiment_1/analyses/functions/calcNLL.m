function nll = calcNLL(responses, p_CW)
    % Handles both single and multiple parameter cases
    %
    % Inputs:
    %   responses: a vector of responses (0 or 1)
    %   p_CW: a vector of pCW values (single param) or matrix (multiple params)
    %
    % Outputs:
    %   nll: the negative log-likelihood (scalar or vector)

    if isvector(p_CW)
        % Single parameter case (current usage)
        ll = responses .* log(p_CW) + (1 - responses) .* log(1 - p_CW);
        nll = -sum(ll, 'omitnan');
    else
        % Multiple parameter case (for vectorization)
        nll = zeros(size(p_CW, 2), 1);
        for i = 1:size(p_CW, 2)
            ll = responses .* log(p_CW(:,i)) + (1 - responses) .* log(1 - p_CW(:,i));
            nll(i) = -sum(ll, 'omitnan');
        end
    end
end 