function r2 = binnedEfronAlongDeltaTheta(delta_theta, y, p_hat, delta_centers, window_width)
% binnedEfronAlongDeltaTheta  Weighted Efron R² on Δθ windows (bin mean y vs bin mean p̂).
%
% For each window (same geometry as mleMuPerDeltaBin), compute n_k, ȳ_k, p̂̄_k, and
%   R² = 1 - sum_k n_k(ȳ_k - p̂̄_k)² / sum_k n_k(ȳ_k - ȳ)²
% with ȳ the global mean response. Summarizes fit of the *trend along Δθ* (like the
% figure), not trial-level 0/1 noise—so values are often much higher than trial Efron.

    y = double(y(:));
    p_hat = double(p_hat(:));
    delta_theta = delta_theta(:);
    n = min([numel(y), numel(p_hat), numel(delta_theta)]);
    y = y(1:n);
    p_hat = p_hat(1:n);
    delta_theta = delta_theta(1:n);

    y_bar_global = mean(y);
    centers = delta_centers(:);

    SS_tot = 0;
    SS_res = 0;

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

        y_b = y(in_win);
        p_b = p_hat(in_win);
        ok = isfinite(y_b) & isfinite(p_b);
        y_b = y_b(ok);
        p_b = p_b(ok);
        if isempty(y_b)
            continue
        end

        n_k = numel(y_b);
        y_bar_k = mean(y_b);
        p_bar_k = mean(p_b);

        SS_tot = SS_tot + n_k * (y_bar_k - y_bar_global)^2;
        SS_res = SS_res + n_k * (y_bar_k - p_bar_k)^2;
    end

    if SS_tot <= eps || ~isfinite(SS_tot) || ~isfinite(SS_res)
        r2 = nan;
    else
        r2 = 1 - SS_res / SS_tot;
    end

end
