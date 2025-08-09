function plotResponseBias(delta, mu, plt_settings, cond, sd_params)

    %% Plot response bias data
    
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