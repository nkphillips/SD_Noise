function plotUnbinnedSdScatterSummaries(fig_dir, summary_table, contrast_lbl, precision_lbl, ps, ci_pct, ci_method)
% plotUnbinnedSdScatterSummaries  Super-Subj-style SD amplitude & FWHM scatter (3 lines × 3 currents).
% Mirrors plotSerialDependence scatter layout (experiment_1/analyses/plotting/plotSerialDependence.m).
%
% Saves PDFs under fig_dir:
%   serial_dependence_amplitude_summary.pdf
%   serial_dependence_fwhm_summary.pdf

    if nargin < 6 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end
    if nargin < 7 || isempty(ci_method)
        ci_method = 'bootstrap';
    end

    pack = packUnbinnedScatterParams(summary_table);
    plot_opts = local_buildPlotOpts(ps);

    % Equal-tailed interval: e.g. [2.5 97.5] -> 95% nominal coverage
    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_lbl = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));

    local_renderOne(fig_dir, pack.amp, pack.amp_lo, pack.amp_hi, ...
        'Amplitude', 'serial_dependence_amplitude_summary.pdf', ...
        sprintf('Serial dependence amplitude (%s)', ci_lbl), ...
        contrast_lbl, precision_lbl, plot_opts, ci_lbl, 'general');

    local_renderOne(fig_dir, pack.fwhm, pack.fwhm_lo, pack.fwhm_hi, ...
        'Width FWHM (deg)', 'serial_dependence_fwhm_summary.pdf', ...
        sprintf('Serial dependence width (FWHM; %s)', ci_lbl), ...
        contrast_lbl, precision_lbl, plot_opts, ci_lbl, 'general');

end

function plot_opts = local_buildPlotOpts(ps)
    plot_opts.line_width = ps.line_width;
    plot_opts.colors = ps.colors;
    plot_opts.figure_color = ps.figure_color;
    plot_opts.tick_length = ps.tick_length;
    plot_opts.marker_size_scatter = 50;
    if isfield(ps, 'marker_size_scatter')
        plot_opts.marker_size_scatter = ps.marker_size_scatter;
    end
    plot_opts.axes_label_font_size = 14;
    plot_opts.axes_tick_font_size = 13;
    try
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    catch
    end
end

function local_renderOne(fig_dir, param_data, ci_lo, ci_hi, ylabel_str, fname, sg_title, ...
    contrast_lbl, precision_lbl, plot_opts, ci_label, ylim_mode)

    has_ci = ~isempty(ci_lo) && ~isempty(ci_hi) && isequal(size(ci_lo), size(param_data)) ...
        && isequal(size(ci_hi), size(param_data)) && any(isfinite(ci_lo(:)) | isfinite(ci_hi(:)));
    y_lims = local_commonYLimits(param_data, ci_lo, ci_hi, has_ci, ylim_mode);

    num_conds = size(param_data, 3);
    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 5]);

    for cond = 1:num_conds

        cond_data = param_data(:, :, cond);

        if cond == 1
            base_color = plot_opts.colors.blue;
            legend_vals = contrast_lbl;
        else
            base_color = plot_opts.colors.green;
            legend_vals = precision_lbl;
        end

        marker_colors = repmat(base_color, size(param_data, 1), 1) .* [1 0.70 0.25]';

        x = 1:3;
        x_dodge = linspace(-0.08, 0.08, size(param_data, 1));
        y = fliplr(cond_data)';

        if has_ci
            ci_lo_cond = ci_lo(:, :, cond);
            ci_hi_cond = ci_hi(:, :, cond);
            lo_y = fliplr(ci_lo_cond)';
            hi_y = fliplr(ci_hi_cond)';
        end

        ax = subplot(1, num_conds, cond);
        hold(ax, 'on');

        for i = 1:size(y, 2)
            xi = x + x_dodge(i);
            plot(ax, xi, y(:, i), '-', 'Color', marker_colors(i, :), 'LineWidth', plot_opts.line_width, 'HandleVisibility', 'off');
            scatter(ax, xi, y(:, i), plot_opts.marker_size_scatter, 'MarkerFaceColor', marker_colors(i, :), ...
                'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.75, 'LineWidth', plot_opts.line_width);
            if has_ci
                yneg = max(0, y(:, i) - lo_y(:, i));
                ypos = max(0, hi_y(:, i) - y(:, i));
                errorbar(ax, xi, y(:, i), yneg, ypos, 'Color', marker_colors(i, :) * 0.8, ...
                    'CapSize', 0, 'LineStyle', 'none', 'LineWidth', plot_opts.line_width, 'HandleVisibility', 'off');
            end
        end

        if cond == 1
            title(ax, 'Contrast', 'FontSize', plot_opts.axes_tick_font_size);
        else
            title(ax, 'Precision', 'FontSize', plot_opts.axes_tick_font_size);
        end

        if cond == 1
            set(ax, 'XTick', 1:3, 'XTickLabel', fliplr(contrast_lbl));
        else
            set(ax, 'XTick', 1:3, 'XTickLabel', fliplr(precision_lbl));
        end

        xlabel(ax, 'Current level', 'FontSize', plot_opts.axes_label_font_size);
        ylabel(ax, ylabel_str, 'FontSize', plot_opts.axes_label_font_size);
        if has_ci && ~isempty(ci_label)
            % Top-left: legend(...,'Location','best') ignores text; bottom-left overlapped CI.
            text(ax, 0.02, 0.98, ci_label, 'Units', 'normalized', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'FontSize', 7, 'Color', plot_opts.colors.gray);
        end

        xlim(ax, [0.75, 3.25]);
        ylim(ax, y_lims);

        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
            'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);

        legend(ax, legend_vals, 'Location', 'best', 'Interpreter', 'none');
    end

    sgtitle(fig, sg_title, 'FontSize', plot_opts.axes_label_font_size + 1);

    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);

end

function y_lims = local_commonYLimits(param_data, ci_lo, ci_hi, has_ci, ylim_mode)
    y_limit_vals = param_data(:);
    if has_ci
        y_limit_vals = [y_limit_vals; ci_lo(:); ci_hi(:)];
    end
    y_limit_vals = y_limit_vals(isfinite(y_limit_vals));

    if isempty(y_limit_vals)
        y_lims = [-1, 1];
        return
    end

    y_min = min(y_limit_vals);
    y_max = max(y_limit_vals);
    if y_min == y_max
        y_pad = max(abs(y_min), 1) * 0.1;
    else
        y_pad = 0.08 * (y_max - y_min);
    end

    if strcmp(ylim_mode, 'nonneg')
        y_lims = [max(0, y_min - y_pad), y_max + y_pad];
    else
        y_lims = [y_min - y_pad, y_max + y_pad];
    end
end

function s = local_ciLabel(ci_method)
    ci_method = lower(char(ci_method));
    switch ci_method
        case 'bca'
            s = 'BCa';
        case 'percentile'
            s = 'percentile';
        otherwise
            s = char(ci_method);
    end
end
