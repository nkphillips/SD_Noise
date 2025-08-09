function plotResponseBias(delta, mu, plt_settings, cond, mu_lo, mu_hi)

    %% Plot response bias data
    
    % Optionally plot shaded CIs if provided
    if nargin >= 6 && ~isempty(mu_lo) && ~isempty(mu_hi)
        % Plot shaded CIs only where we have finite bounds
        delta = delta(:)';
        mu = mu(:)';
        mu_lo = mu_lo(:)';
        mu_hi = mu_hi(:)';
        idx = isfinite(mu) & isfinite(mu_lo) & isfinite(mu_hi);
        if any(idx)
            err_upper = mu_hi(idx) - mu(idx);
            err_lower = mu(idx) - mu_lo(idx);
            err = [err_upper; err_lower];
            if cond == 1
                lineProps = {'-','Color', plt_settings.colors.blue, 'LineWidth', plt_settings.line_width};
                lineColor = plt_settings.colors.blue;
            else
                lineProps = {'-','Color', plt_settings.colors.green, 'LineWidth', plt_settings.line_width};
                lineColor = plt_settings.colors.green;
            end
            shadedErrorBar(delta(idx), mu(idx), err, 'lineProps', lineProps, 'transparent', true, 'patchSaturation', 0.5);
            hold on;
            % Also overlay the full line for context
            plot(delta, mu, 'LineWidth', plt_settings.line_width, 'Color', lineColor);
            hold on;
            return
        end
        % If no finite CIs, fall through to plain line plot below
    end

    % Use different colors for different conditions
    if cond == 1
        % Contrast condition - blue
        plot(delta, mu, 'LineWidth', plt_settings.line_width, 'Color', plt_settings.colors.blue);
    elseif cond == 2
        % Precision condition - green
        plot(delta, mu, 'LineWidth', plt_settings.line_width, 'Color', plt_settings.colors.green);
    end

    hold on;
end