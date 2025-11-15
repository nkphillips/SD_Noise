%% plotSerialDependence
% Bar plots of serial dependence estimates (amplitude and width)
% Each subplot in the figure is the noise type, bar group is the previous level, while each bar is the current level.

% Inputs:
% p: parameters
% plt_opts: plotting settings
% sd: serial dependence estimates

% Outputs:
% fg: figure handle

function plotSerialDependence(sd_data, param_index, param_name, p, plt_opts, fg, sd_ci_lo, sd_ci_hi)
% plotSerialDependence - Create heatmap of serial dependence parameters
%
% This function creates a heatmap visualization of serial dependence parameters
% across all experimental conditions (3×3×2 structure).
%
% Inputs:
%   sd_data - 4D array [prev_lvl, curr_lvl, cond, param] from sd.all.params_est
%   param_index - which parameter to plot (1=amplitude, 2=width, 3=baseline, 4=sigma)
%   param_name - string name of parameter for title
%   p - parameter struct
%   plt_opts - plotting settings struct
%   sd_ci_lo - optional 4D array [prev, curr, cond, param] lower CI
%   sd_ci_hi - optional 4D array [prev, curr, cond, param] upper CI

% Extract the parameter data
param_data = squeeze(sd_data(:,:,:,param_index)); % [prev_lvl, curr_lvl, cond]

% If plotting width (param_index==2), convert from model parameter w (1/deg)
% to Gaussian FWHM in degrees: FWHM = (2 * sqrt(ln 2)) / w ≈ 1.6651 / w
if param_index == 2
    with_epsilon = max(param_data, eps);
    param_data = (2 * sqrt(log(2))) ./ with_epsilon;
    % Clarify label for width
    param_name = [param_name ' (FWHM, °)'];
end

% Extract CIs if provided (match param and apply FWHM transform for width)
has_ci = (nargin >= 8) && ~isempty(sd_ci_lo) && ~isempty(sd_ci_hi);
if has_ci
    ci_lo = squeeze(sd_ci_lo(:,:,:,param_index));
    ci_hi = squeeze(sd_ci_hi(:,:,:,param_index));
    if param_index == 2
        % Transform width CIs to FWHM: FWHM = c / w, with monotone decreasing
        c = 2 * sqrt(log(2));
        with_eps_lo = max(ci_lo, eps);
        with_eps_hi = max(ci_hi, eps);
        % Since c/w is decreasing, lo_FWHM = c/hi_w, hi_FWHM = c/lo_w
        lo_fwhm = c ./ with_eps_hi;
        hi_fwhm = c ./ with_eps_lo;
        ci_lo = lo_fwhm;
        ci_hi = hi_fwhm;
    end
end

% Create figure with subplots for each condition
num_conds = size(param_data, 3);

% Calculate global min and max for consistent colorbar range when needed
if plt_opts.sd_colorbar_global
    global_min = min(param_data(:));
    global_max = max(param_data(:));
    if global_min == global_max
        % Avoid zero-range color limits
        pad = max(abs(global_min), 1) * eps;
        global_min = global_min - pad;
        global_max = global_max + pad;
    end
end

%% Create heatmap with chosen colorbar scaling

for cond = 1:num_conds

    subplot(1, num_conds, cond);

    % Get data for current condition
    cond_data = param_data(:,:,cond);

    if plt_opts.sd_colorbar_global
        imagesc(cond_data, [global_min, global_max]);
    else
        local_min = min(cond_data(:));
        local_max = max(cond_data(:));
        if local_min == local_max
            pad = max(abs(local_min), 1) * eps;
            local_min = local_min - pad;
            local_max = local_max + pad;
        end
        imagesc(cond_data, [local_min, local_max]);
    end

    % Condition-specific colormap: low values = dark hue, high values = light of same hue
    if cond == 1
        base_color = plt_opts.colors.blue; % contrast
    else
        base_color = plt_opts.colors.green; % precision
    end
    start_color = 0.35 * base_color; % dark tone
    white_color = [1 1 1];
    steps = 256;
    cm = [linspace(start_color(1), white_color(1), steps)', ...
        linspace(start_color(2), white_color(2), steps)', ...
        linspace(start_color(3), white_color(3), steps)'];
    cm = flipud(cm);
    colormap(gca, cm);

    % Add colorbar
    colorbar;

    % Set axis labels and ticks
    xlabel('Current Level');
    ylabel('Previous Level');

    % Use actual contrast/precision values for tick labels
    if cond == 1
        % Contrast condition
        set(gca, 'XTick', 1:3, 'XTickLabel', p.contrast);
        set(gca, 'YTick', 1:3, 'YTickLabel', p.contrast);
    else
        % Precision condition
        set(gca, 'XTick', 1:3, 'XTickLabel', p.precision);
        set(gca, 'YTick', 1:3, 'YTickLabel', p.precision);
    end

    % Add value annotations
    for i = 1:3
        for j = 1:3
            text(j, i, sprintf('%.2f', cond_data(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', plt_opts.colors.gray - 0.1, 'FontWeight', 'bold', ...
                'FontSize', 10);
        end
    end

    % Set title
    if cond == 1
        title('Contrast');
    else
        title('Precision');
    end

    % Set axis properties
    axis square;
    set(gca, 'TickDir', 'out', 'LineWidth', plt_opts.line_width, 'Box', 'off');

end

% Save figure
if plt_opts.save_sup_figures
    fg_name = ['Super Subj SD ' param_name ' Heat Map'];
    saveas(fg, fullfile(plt_opts.sup_figure_path, [fg_name '.' plt_opts.fg_type]));
end

clf(fg);

%% Create scatter plot

for cond = 1:num_conds

    % Get data for current condition
    cond_data = param_data(:,:,cond);

    % Set Colors
    if cond == 1
        base_color = plt_opts.colors.blue; % contrast
        legend_vals = p.contrast;
    else
        base_color = plt_opts.colors.green; % precision
        legend_vals = p.precision;
    end

    marker_colors = repmat(base_color, size(param_data,1),1) .* [1 0.70 0.25]';

    % Get local y min and max
    local_min = min(cond_data(:));
    local_max = max(cond_data(:));
    if local_min == local_max
        pad = max(abs(local_min), 1) * eps;
        local_min = local_min - pad;
        local_max = local_max + pad;
    end

    % Get global y min and max
    global_min = min(param_data(:));
    global_max = max(param_data(:));
    if global_min == global_max
        pad = max(abs(global_min), 1) * eps;
        global_min = global_min - pad;
        global_max = global_max + pad;
    end

    % Round to multiples of 5
    local_min = floor(local_min / 5) * 5;
    local_max = ceil(local_max / 5) * 5;

    global_min = floor(global_min / 5) * 5;
    global_max = ceil(global_max / 5) * 5;

    % Plot data
    subplot(1, num_conds, cond)

    x = 1:3;
    y = fliplr(cond_data)';
    hold on

    % Prepare CI for this condition if available
    if has_ci
        ci_lo_cond = squeeze(ci_lo(:,:,cond));
        ci_hi_cond = squeeze(ci_hi(:,:,cond));
        lo_y = fliplr(ci_lo_cond)';
        hi_y = fliplr(ci_hi_cond)';
    end

    for i = 1:size(y,2)
        plot(x, y(:,i), '-', 'Color', marker_colors(i,:), 'LineWidth', plt_opts.line_width, 'HandleVisibility', 'off')
        scatter(x, y(:,i), 50, 'MarkerFaceColor', marker_colors(i,:), 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.75, 'LineWidth', plt_opts.line_width)
        if has_ci && (param_index == 1 || param_index == 2)
            % Asymmetric error bars per point
            yneg = max(0, y(:,i) - lo_y(:,i));
            ypos = max(0, hi_y(:,i) - y(:,i));
            errorbar(x, y(:,i), yneg, ypos, 'Color', marker_colors(i,:)*0.8, 'CapSize', 0, 'LineStyle', 'none', 'LineWidth', plt_opts.line_width, 'HandleVisibility', 'off');
        end
    end

    % Format figure
    if cond == 1
        title('Contrast');
    else
        title('Precision');
    end

    if cond == 1
        % Contrast condition
        set(gca, 'XTick', 1:3, 'XTickLabel', fliplr(p.contrast));
    else
        % Precision condition
        set(gca, 'XTick', 1:3, 'XTickLabel', fliplr(p.precision));
    end

    xlabel('Current level')
    ylabel(param_name)

    % ylim([local_min local_max])
    % ylim([global_min global_max])
    if strcmp(param_name, 'Amplitude'), ylim([0 10]); end

    axis square;
    set(gca, 'TickDir', 'out', 'LineWidth', plt_opts.line_width, 'Box', 'off');


    legend(legend_vals)

end

% Save figure
if plt_opts.save_sup_figures
    fg_name = ['Super Subj SD ' param_name ' Scatter'];
    saveas(fg, fullfile(plt_opts.sup_figure_path, [fg_name '.' plt_opts.fg_type]));
end

clf(fg);

end