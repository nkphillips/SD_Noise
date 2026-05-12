function plotUnbinnedSdScatterByNBack(fig_dir, results, n_back_list, contrast_lbl, precision_lbl, ps, ci_pct, ci_method)
% plotUnbinnedSdScatterByNBack  Scatter summaries split by manipulation across n-back.
%
% For each metric (amplitude, FWHM) and manipulation (contrast, precision),
% renders one figure with one subplot per n-back. Axes are shared within each
% figure, making n-back changes directly comparable.

    if nargin < 7 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end
    if nargin < 8 || isempty(ci_method)
        ci_method = local_resultsCIMethod(results);
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    plot_opts = local_buildPlotOpts(ps);
    n_back_list = n_back_list(:)';
    n_nb = numel(n_back_list);
    num_levels = 3;
    num_conds = 2;

    amp = nan(num_levels, num_levels, num_conds, n_nb);
    amp_lo = amp;
    amp_hi = amp;
    fwhm = amp;
    fwhm_lo = amp;
    fwhm_hi = amp;

    for i_nb = 1:n_nb
        res_i = local_getResult(results, n_back_list(i_nb), i_nb);
        if isempty(res_i) || ~isfield(res_i, 'summary_table') || isempty(res_i.summary_table)
            continue
        end
        pack = packUnbinnedScatterParams(res_i.summary_table);
        amp(:, :, :, i_nb) = pack.amp;
        amp_lo(:, :, :, i_nb) = pack.amp_lo;
        amp_hi(:, :, :, i_nb) = pack.amp_hi;
        fwhm(:, :, :, i_nb) = pack.fwhm;
        fwhm_lo(:, :, :, i_nb) = pack.fwhm_lo;
        fwhm_hi(:, :, :, i_nb) = pack.fwhm_hi;
    end

    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));

    local_renderMetric(fig_dir, amp, amp_lo, amp_hi, n_back_list, contrast_lbl, precision_lbl, ...
        plot_opts, ci_label, 'Amplitude A (deg)', ...
        'amplitude', 'Unbinned MLE: SD amplitude across n-back', false);

    local_renderMetric(fig_dir, fwhm, fwhm_lo, fwhm_hi, n_back_list, contrast_lbl, precision_lbl, ...
        plot_opts, ci_label, 'Width FWHM (deg)', ...
        'width_fwhm', 'Unbinned MLE: SD width (FWHM) across n-back', true);
end

function local_renderMetric(fig_dir, vals, ci_lo, ci_hi, n_back_list, contrast_lbl, precision_lbl, ...
    plot_opts, ci_label, ylabel_str, metric_file, title_base, nonneg_y)

    for cond = 1:2
        if cond == 1
            labels = contrast_lbl;
            cond_name = 'Contrast';
            cond_file = 'contrast';
            base_color = plot_opts.colors.blue;
        else
            labels = precision_lbl;
            cond_name = 'Precision';
            cond_file = 'precision';
            base_color = plot_opts.colors.green;
        end

        y_lims = local_commonYLimits(vals(:, :, cond, :), ci_lo(:, :, cond, :), ...
            ci_hi(:, :, cond, :), nonneg_y);
        fig_w = max(6, 3.6 * numel(n_back_list));
        fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
            'Units', 'inches', 'Position', [1 1 fig_w 4.5], ...
            'PaperUnits', 'inches', 'PaperSize', [fig_w 4.5], ...
            'PaperPositionMode', 'manual', 'PaperPosition', [0 0 fig_w 4.5]);
        tl = tiledlayout(1, numel(n_back_list), 'Padding', 'compact', 'TileSpacing', 'compact');

        marker_colors = repmat(base_color, 3, 1) .* [1 0.70 0.25]';
        x = 1:3;
        x_dodge = linspace(-0.08, 0.08, 3);

        for i_nb = 1:numel(n_back_list)
            ax = nexttile(tl);
            hold(ax, 'on');

            cond_data = vals(:, :, cond, i_nb);
            y = fliplr(cond_data)';
            lo_y = fliplr(ci_lo(:, :, cond, i_nb))';
            hi_y = fliplr(ci_hi(:, :, cond, i_nb))';

            for prev = 1:size(y, 2)
                xi = x + x_dodge(prev);
                plot(ax, xi, y(:, prev), '-', 'Color', marker_colors(prev, :), ...
                    'LineWidth', plot_opts.line_width, 'HandleVisibility', 'off');
                scatter(ax, xi, y(:, prev), plot_opts.marker_size_scatter, ...
                    'MarkerFaceColor', marker_colors(prev, :), 'MarkerEdgeColor', [1 1 1], ...
                    'MarkerFaceAlpha', 0.75, 'LineWidth', plot_opts.line_width, ...
                    'DisplayName', labels{prev});

                yneg = max(0, y(:, prev) - lo_y(:, prev));
                ypos = max(0, hi_y(:, prev) - y(:, prev));
                errorbar(ax, xi, y(:, prev), yneg, ypos, ...
                    'Color', marker_colors(prev, :) * 0.8, 'CapSize', 0, ...
                    'LineStyle', 'none', 'LineWidth', plot_opts.line_width, ...
                    'HandleVisibility', 'off');
            end

            if ~nonneg_y
                yline(ax, 0, 'k-', 'HandleVisibility', 'off');
            end
            xlim(ax, [0.75, 3.25]);
            ylim(ax, y_lims);
            set(ax, 'XTick', 1:3, 'XTickLabel', fliplr(labels));
            xlabel(ax, 'Current level', 'FontSize', plot_opts.axes_label_font_size);
            if i_nb == 1
                ylabel(ax, ylabel_str, 'FontSize', plot_opts.axes_label_font_size, 'Interpreter', 'tex');
            end
            title(ax, sprintf('%d-back', n_back_list(i_nb)), ...
                'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
            axis(ax, 'square');
            box(ax, 'off');
            set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
                'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);

            if i_nb == numel(n_back_list)
                legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
            end
        end

        title(tl, sprintf('%s: %s (%s)', title_base, cond_name, ci_label), ...
            'FontSize', plot_opts.axes_label_font_size + 1);
        out_pdf = fullfile(fig_dir, sprintf('unbinned_mle_super_sd_%s_scatter_by_nback_%s.pdf', ...
            metric_file, cond_file));
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
        close(fig);
    end
end

function y_lims = local_commonYLimits(vals, ci_lo, ci_hi, nonneg_y)
    y_limit_vals = [vals(:); ci_lo(:); ci_hi(:)];
    y_limit_vals = y_limit_vals(isfinite(y_limit_vals));

    if isempty(y_limit_vals)
        y_lims = [-1, 1];
        if nonneg_y
            y_lims = [0, 1];
        end
        return
    end

    y_min = min(y_limit_vals);
    y_max = max(y_limit_vals);
    if y_min == y_max
        y_pad = max(abs(y_min), 1) * 0.1;
    else
        y_pad = 0.08 * (y_max - y_min);
    end

    if nonneg_y
        y_lims = [max(0, y_min - y_pad), y_max + y_pad];
    else
        y_lims = [y_min - y_pad, y_max + y_pad];
    end
end

function res_i = local_getResult(results, n_back, i_nb)
    if isstruct(results) && isfield(results, sprintf('n%d', n_back))
        res_i = results.(sprintf('n%d', n_back));
    elseif isstruct(results) && isfield(results, 'summary_table') && i_nb == 1
        res_i = results;
    else
        res_i = [];
    end
end

function ci_method = local_resultsCIMethod(results)
    ci_method = 'bootstrap';
    if isstruct(results) && isfield(results, 'ci_method') && ~isempty(results.ci_method)
        ci_method = results.ci_method;
        return
    end
    if isstruct(results)
        fields = fieldnames(results);
        for i = 1:numel(fields)
            val = results.(fields{i});
            if isstruct(val) && isfield(val, 'ci_method') && ~isempty(val.ci_method)
                ci_method = val.ci_method;
                return
            end
        end
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
