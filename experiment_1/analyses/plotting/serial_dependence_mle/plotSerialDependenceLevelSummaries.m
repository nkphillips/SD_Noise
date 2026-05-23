function out = plotSerialDependenceLevelSummaries(fig_dir, summary_table, contrast_lbl, precision_lbl, ps, ci_pct, ci_method, n_back)
% plotSerialDependenceLevelSummaries  Single-n-back 1x2 level-summary figures.

    if nargin < 6 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end
    if nargin < 7 || isempty(ci_method)
        ci_method = 'bootstrap';
    end
    if nargin < 8 || isempty(n_back)
        n_back = NaN;
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    pack = packSerialDependenceScatterParams(summary_table);
    plot_opts = local_buildPlotOpts(ps);
    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));
    nback_tag = local_nBackTag(n_back);

    out = struct();
    out.amplitude_by_past_pdf = local_renderByAxis(fig_dir, pack.amp, pack.amp_lo, pack.amp_hi, ...
        'Amplitude (deg)', sprintf('serial_dependence_amplitude_by_%s_by_past_level.pdf', nback_tag), ...
        sprintf('Serial dependence amplitude, %s, by past level (%s)', local_nBackLabel(n_back), ci_label), ...
        contrast_lbl, precision_lbl, plot_opts, 'past');
    out.amplitude_by_current_pdf = local_renderByAxis(fig_dir, pack.amp, pack.amp_lo, pack.amp_hi, ...
        'Amplitude (deg)', sprintf('serial_dependence_amplitude_by_%s_by_current_level.pdf', nback_tag), ...
        sprintf('Serial dependence amplitude, %s, by current level (%s)', local_nBackLabel(n_back), ci_label), ...
        contrast_lbl, precision_lbl, plot_opts, 'current');
    out.fwhm_by_past_pdf = local_renderByAxis(fig_dir, pack.fwhm, pack.fwhm_lo, pack.fwhm_hi, ...
        'FWHM (deg)', sprintf('serial_dependence_fwhm_by_%s_by_past_level.pdf', nback_tag), ...
        sprintf('Serial dependence FWHM, %s, by past level (%s)', local_nBackLabel(n_back), ci_label), ...
        contrast_lbl, precision_lbl, plot_opts, 'past');
    out.fwhm_by_current_pdf = local_renderByAxis(fig_dir, pack.fwhm, pack.fwhm_lo, pack.fwhm_hi, ...
        'FWHM (deg)', sprintf('serial_dependence_fwhm_by_%s_by_current_level.pdf', nback_tag), ...
        sprintf('Serial dependence FWHM, %s, by current level (%s)', local_nBackLabel(n_back), ci_label), ...
        contrast_lbl, precision_lbl, plot_opts, 'current');
end

function out_pdf = local_renderByAxis(fig_dir, param_data, ci_lo, ci_hi, ylabel_str, fname, sg_title, ...
    contrast_lbl, precision_lbl, plot_opts, x_axis)

    has_ci = ~isempty(ci_lo) && ~isempty(ci_hi) && isequal(size(ci_lo), size(param_data)) ...
        && isequal(size(ci_hi), size(param_data)) && any(isfinite(ci_lo(:)) | isfinite(ci_hi(:)));
    y_lims = local_commonYLimits(param_data, ci_lo, ci_hi, has_ci);

    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 5]);

    cond_order = [2 1];   % precision first, contrast second
    for panel = 1:numel(cond_order)
        cond = cond_order(panel);
        if cond == 1
            base_color = plot_opts.colors.blue;
            line_labels = contrast_lbl;
            title_str = 'Contrast';
        else
            base_color = plot_opts.colors.green;
            line_labels = precision_lbl;
            title_str = 'Precision';
        end

        marker_colors = repmat(base_color, size(param_data, 1), 1) .* [1 0.70 0.25]';
        ax = subplot(1, 2, panel);
        hold(ax, 'on');
        x = 1:3;
        x_dodge = linspace(-0.08, 0.08, size(param_data, 1));

        for line_level = 1:3
            if strcmp(x_axis, 'past')
                y = squeeze(param_data(:, line_level, cond));
                lo_y = squeeze(ci_lo(:, line_level, cond));
                hi_y = squeeze(ci_hi(:, line_level, cond));
                display_name = sprintf('Current %s', line_labels{line_level});
            else
                y = squeeze(param_data(line_level, :, cond));
                lo_y = squeeze(ci_lo(line_level, :, cond));
                hi_y = squeeze(ci_hi(line_level, :, cond));
                display_name = sprintf('Past %s', line_labels{line_level});
            end
            xi = x + x_dodge(line_level);
            plot(ax, xi, y, '-', 'Color', marker_colors(line_level, :), ...
                'LineWidth', plot_opts.line_width, 'HandleVisibility', 'off');
            scatter(ax, xi, y, plot_opts.marker_size_scatter, ...
                'MarkerFaceColor', marker_colors(line_level, :), ...
                'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.75, ...
                'LineWidth', plot_opts.line_width, 'DisplayName', display_name);
            if has_ci
                yneg = max(0, y - lo_y);
                ypos = max(0, hi_y - y);
                errorbar(ax, xi, y, yneg, ypos, 'Color', marker_colors(line_level, :) * 0.8, ...
                    'CapSize', 0, 'LineStyle', 'none', 'LineWidth', plot_opts.line_width, ...
                    'HandleVisibility', 'off');
            end
        end

        title(ax, title_str, 'FontSize', plot_opts.axes_tick_font_size);
        if strcmp(x_axis, 'past')
            xlabel(ax, 'Past level', 'FontSize', plot_opts.axes_label_font_size);
        else
            xlabel(ax, 'Current level', 'FontSize', plot_opts.axes_label_font_size);
        end
        ylabel(ax, ylabel_str, 'FontSize', plot_opts.axes_label_font_size);
        set(ax, 'XTick', 1:3, 'XTickLabel', local_levelLabels(cond, contrast_lbl, precision_lbl));
        xlim(ax, [0.75, 3.25]);
        ylim(ax, y_lims);
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
            'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
        legend(ax, 'Location', 'best', 'Interpreter', 'none');
    end

    sgtitle(fig, sg_title, 'FontSize', plot_opts.axes_label_font_size + 1, 'Interpreter', 'none');
    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function labels = local_levelLabels(cond, contrast_lbl, precision_lbl)
    if cond == 1
        labels = contrast_lbl;
    else
        labels = precision_lbl;
    end
    if isstring(labels)
        labels = cellstr(labels);
    end
    labels = fliplr(labels);
end

function y_lims = local_commonYLimits(param_data, ci_lo, ci_hi, has_ci)
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
    y_lims = [y_min - y_pad, y_max + y_pad];
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

function tag = local_nBackTag(n_back)
    if isfinite(n_back)
        tag = sprintf('%dback', n_back);
    else
        tag = 'nback';
    end
end

function label = local_nBackLabel(n_back)
    if isfinite(n_back)
        label = sprintf('%d-back', n_back);
    else
        label = 'n-back';
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
