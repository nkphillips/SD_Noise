function p_CW = calc_pCW(x, mu, sigma, guess_rate)
    % calc_pCW  Predicted P(CW) under a Gaussian psychometric with constant lapse rate.
    %
    %   P(CW) = (1 - guess_rate) * normcdf(x, mu, sigma) + 0.5 * guess_rate
    %
    % Three call patterns are supported:
    %   1) Scalar (mu, sigma): one psychometric for all probe offsets x.
    %   2) Trial-aligned vector mu of length N == numel(x), scalar (or length-N) sigma:
    %      computes p_CW element-wise. This is what calcSerialDependenceFit needs in NLL
    %      mode where mu_t = DoG(delta_theta_t) + baseline.
    %   3) Multiple parameter pairs (length-K mu, length-K sigma) with K ~= numel(x):
    %      returns an N x K matrix evaluating each pair at every x (used for grid search
    %      vectorization, though not currently called this way).
    %
    % Returns p_CW clipped to >= 1e-10 to avoid log(0) downstream.

    x = x(:);
    nx = numel(x);

    if isscalar(mu) && isscalar(sigma)
        % Case 1: scalar parameters
        p_CW = (1 - guess_rate) * normcdf(x, mu, sigma) + 0.5 * guess_rate;
    elseif numel(mu) == nx && (isscalar(sigma) || numel(sigma) == nx)
        % Case 2: trial-aligned (per-trial mu, scalar or per-trial sigma)
        mu = mu(:);
        if ~isscalar(sigma), sigma = sigma(:); end
        p_CW = (1 - guess_rate) * normcdf(x, mu, sigma) + 0.5 * guess_rate;
    else
        % Case 3: multiple (mu_i, sigma_i) parameter pairs
        n_pairs = numel(mu);
        p_CW = zeros(nx, n_pairs);
        for i = 1:n_pairs
            p_CW(:, i) = (1 - guess_rate) * normcdf(x, mu(i), sigma(i)) + 0.5 * guess_rate;
        end
    end

    % Avoid log(0)
    p_CW = max(p_CW, 1e-10);

end
