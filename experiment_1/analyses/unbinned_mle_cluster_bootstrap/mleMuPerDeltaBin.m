function [delta_centers, mu_deg] = mleMuPerDeltaBin(delta_theta, x_probe, response, sigma, delta_centers, window_width, guess_rate)
% mleMuPerDeltaBin  Per Δθ window: scalar MLE of psychometric location μ with fixed σ.
%
%   P(CW) = (1 - lambda) * normcdf(x_probe, mu, sigma) + 0.5 * lambda     (S&S 2022 lapse model)
%
% Same likelihood as fitConditionMLE; σ fixed (typically pooled σ̂) so each bin has
% one free μ comparable to μ(Δθ) = DoG(Δθ)+β on the figure y-axis.
% Window logic matches empiricalBiasProxyBinned (wrap at ±90°).

    if nargin < 7 || isempty(guess_rate)
        guess_rate = 0.25;
    end

    delta_theta = delta_theta(:);
    x_probe = x_probe(:);
    response = double(response(:));

    centers = delta_centers(:);
    mu_deg = nan(size(centers));

    for i = 1:numel(centers)
        c = centers(i);
        left_edge = c - window_width / 2;
        right_edge = c + window_width / 2;

        if left_edge < -90
            in_win = delta_theta >= (left_edge + 180) | delta_theta <= right_edge;
        elseif right_edge > 90
            in_win = delta_theta <= (right_edge - 180) | delta_theta >= left_edge;
        else
            in_win = delta_theta >= left_edge & delta_theta <= right_edge;
        end

        if ~any(in_win)
            continue
        end

        xp = x_probe(in_win);
        rv = response(in_win);
        valid = rv == 0 | rv == 1;
        xp = xp(valid);
        rv = rv(valid);

        if isempty(xp)
            continue
        end

        mu_deg(i) = mleMuScalarBin(xp, rv, sigma, guess_rate);
    end

    delta_centers = centers;

end

function mu_hat = mleMuScalarBin(x_probe, response, sigma, guess_rate)
% One-dimensional MLE for μ; σ fixed; 25% guess rate.

    x_probe = double(x_probe(:));
    response = double(response(:));
    n = min(numel(x_probe), numel(response));
    x_probe = x_probe(1:n);
    response = response(1:n);

    if numel(unique(response)) < 2
        mu_hat = nan;
        return
    end

    eps_log = eps;
    nll = @(mu) -sum( ...
        response .* log(((1 - guess_rate) .* normcdf(x_probe, mu, sigma) + 0.5 .* guess_rate) + eps_log) + ...
        (1 - response) .* log(1 - ((1 - guess_rate) .* normcdf(x_probe, mu, sigma) + 0.5 .* guess_rate) + eps_log), ...
        'omitnan');

    try
        mu_hat = fminbnd(nll, -45, 45);
        if ~isfinite(mu_hat)
            mu_hat = nan;
        end
    catch
        mu_hat = nan;
    end

end
