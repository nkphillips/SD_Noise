function [delta_centers, bias_deg] = empiricalBiasProxyBinned(delta_theta, x_probe, response, sigma, delta_centers, window_width)
% empiricalBiasProxyBinned  Coarse empirical "bias" vs Δθ for overlay plots (visualization only).
%
% For each Δθ window, p_CW = mean(response). Map to a shift in degrees using the same σ as the
% psychometric link: bias ≈ σ * (norminv(p_CW) - norminv(0.5)), clamped to finite values.
% Matches ordering of delta_theta_centers with wrap handling like makeDeltaThetaWindows.

    delta_theta = delta_theta(:);
    x_probe = x_probe(:);
    response = double(response(:));

    centers = delta_centers(:);
    bias_deg = nan(size(centers));

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

        r = response(in_win);
        r = r(r == 0 | r == 1);
        if isempty(r)
            continue
        end

        p_cw = mean(r);
        p_cw = min(max(p_cw, eps), 1 - eps);
        bias_deg(i) = sigma * (norminv(p_cw) - norminv(0.5));
    end

    delta_centers = centers;

end
