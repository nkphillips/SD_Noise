function plotSerialDependencePooledSubjectPointSummaries(fig_dir, summary_table, subject_cell_fits, contrast_lbl, precision_lbl, ps)
% plotSerialDependencePooledSubjectPointSummaries  Pooled estimates with faint subject points.

    if isempty(subject_cell_fits) || height(subject_cell_fits) == 0
        return
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    pack = packSerialDependenceScatterParams(summary_table);
    subj_pack = local_packSubjectFits(subject_cell_fits);
    plot_opts = local_buildPlotOpts(ps);

    local_renderOne(fig_dir, pack.amp, subj_pack.amp, subj_pack.at_bound, ...
        'Amplitude', 'serial_dependence_amplitude_subject_points.pdf', ...
        'Serial dependence amplitude with subject-level points', ...
        contrast_lbl, precision_lbl, plot_opts);

    local_renderOne(fig_dir, pack.fwhm, subj_pack.fwhm, subj_pack.at_bound, ...
        'Width FWHM (deg)', 'serial_dependence_fwhm_subject_points.pdf', ...
        'Serial dependence FWHM with subject-level points', ...
        contrast_lbl, precision_lbl, plot_opts);
end

function subj_pack = local_packSubjectFits(T)
    subj_list = unique(T.subject_id, 'stable');
    n_subj = numel(subj_list);
    subj_pack.amp = nan(n_subj, 3, 3, 2);
    subj_pack.fwhm = nan(n_subj, 3, 3, 2);
    subj_pack.at_bound = false(n_subj, 3, 3, 2);

    for r = 1:height(T)
        s = find(subj_list == T.subject_id(r), 1);
        if isempty(s)
            continue
        end
        prev = T.cond_prev(r);
        curr = T.cond_curr(r);
        if strcmpi(char(T.manipulation(r)), 'contrast')
            m = 1;
        else
            m = 2;
        end
        subj_pack.amp(s, prev, curr, m) = T.A(r);
        subj_pack.fwhm(s, prev, curr, m) = T.fwhm(r);
        if ismember('at_bound_any', T.Properties.VariableNames)
            subj_pack.at_bound(s, prev, curr, m) = logical(T.at_bound_any(r));
        end
    end
end

function plot_opts = local_buildPlotOpts(ps)
    plot_opts.line_width = ps.line_width;
    plot_opts.colors = ps.colors;
    plot_opts.figure_color = ps.figure_color;
    plot_opts.tick_length = ps.tick_length;
    plot_opts.marker_size_scatter = 55;
    if isfield(ps, 'marker_size_scatter')
        plot_opts.marker_size_scatter = ps.marker_size_scatter;
    end
    plot_opts.axes_label_font_size = 14;
    plot_opts.axes_tick_font_size = 13;
    if isfield(ps, 'axes_label_font_size')
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
    end
    if isfield(ps, 'axes_tick_font_size')
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    end
end

function local_renderOne(fig_dir, pooled_data, subj_data, subj_at_bound, ylabel_str, fname, sg_title, ...
    contrast_lbl, precision_lbl, plot_opts)

    y_lims = local_commonYLimits(pooled_data, subj_data);
    num_conds = size(pooled_data, 3);

    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 5]);

    cond_order = [2 1];   % precision first, contrast second
    for panel = 1:num_conds
        cond = cond_order(panel);
        if cond == 1
            base_color = plot_opts.colors.blue;
            legend_vals = contrast_lbl;
            title_str = 'Contrast';
        else
            base_color = plot_opts.colors.green;
            legend_vals = precision_lbl;
            title_str = 'Precision';
        end

        marker_colors = repmat(base_color, size(pooled_data, 1), 1) .* [1 0.70 0.25]';
        x = 1:3;
        x_dodge = linspace(-0.08, 0.08, size(pooled_data, 1));
        pooled_y = fliplr(pooled_data(:, :, cond))';

        ax = subplot(1, num_conds, panel);
        hold(ax, 'on');

        for s = 1:size(subj_data, 1)
            subj_mat = squeeze(subj_data(s, :, :, cond));
            subj_bound = squeeze(subj_at_bound(s, :, :, cond));
            subj_y = fliplr(subj_mat)';
            bound_y = fliplr(subj_bound)';
            subj_jitter = local_subjectJitter(s, size(subj_data, 1));
            for i_prev = 1:size(pooled_y, 2)
                xi = x + x_dodge(i_prev) + subj_jitter;
                good = isfinite(subj_y(:, i_prev)) & ~bound_y(:, i_prev);
                bounded = isfinite(subj_y(:, i_prev)) & bound_y(:, i_prev);
                if any(good)
                    scatter(ax, xi(good), subj_y(good, i_prev), 18, marker_colors(i_prev, :), ...
                        'filled', 'MarkerFaceAlpha', 0.16, 'MarkerEdgeAlpha', 0.08, ...
                        'HandleVisibility', 'off');
                end
                if any(bounded)
                    scatter(ax, xi(bounded), subj_y(bounded, i_prev), 18, marker_colors(i_prev, :), ...
                        'MarkerFaceColor', 'none', 'MarkerEdgeAlpha', 0.20, ...
                        'HandleVisibility', 'off');
                end
            end
        end

        for i_prev = 1:size(pooled_y, 2)
            xi = x + x_dodge(i_prev);
            plot(ax, xi, pooled_y(:, i_prev), '-', 'Color', marker_colors(i_prev, :) * 0.75, ...
                'LineWidth', plot_opts.line_width + 0.5, 'HandleVisibility', 'off');
            scatter(ax, xi, pooled_y(:, i_prev), plot_opts.marker_size_scatter, ...
                'MarkerFaceColor', marker_colors(i_prev, :), 'MarkerEdgeColor', [1 1 1], ...
                'MarkerFaceAlpha', 0.95, 'LineWidth', plot_opts.line_width, ...
                'DisplayName', legend_vals{i_prev});
        end

        title(ax, title_str, 'FontSize', plot_opts.axes_tick_font_size);
        if cond == 1
            set(ax, 'XTick', 1:3, 'XTickLabel', fliplr(contrast_lbl));
        else
            set(ax, 'XTick', 1:3, 'XTickLabel', fliplr(precision_lbl));
        end
        xlabel(ax, 'Current level', 'FontSize', plot_opts.axes_label_font_size);
        ylabel(ax, ylabel_str, 'FontSize', plot_opts.axes_label_font_size);
        text(ax, 0.02, 0.98, 'Faint points = subject fits; bold markers = pooled fit', ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Color', plot_opts.colors.gray);
        xlim(ax, [0.75, 3.25]);
        ylim(ax, y_lims);
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
            'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
        legend(ax, legend_vals, 'Location', 'best', 'Interpreter', 'none');
    end

    sgtitle(fig, sg_title, 'FontSize', plot_opts.axes_label_font_size + 1);
    exportgraphics(fig, fullfile(fig_dir, fname), 'ContentType', 'vector');
    close(fig);
end

function jitter = local_subjectJitter(s, n_subj)
    if n_subj <= 1
        jitter = 0;
    else
        jitter = -0.025 + 0.05 * (s - 1) / (n_subj - 1);
    end
end

function y_lims = local_commonYLimits(pooled_data, subj_data)
    vals = [pooled_data(:); subj_data(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_lims = [-1, 1];
        return
    end
    y_min = min(vals);
    y_max = max(vals);
    if y_min == y_max
        y_pad = max(abs(y_min), 1) * 0.1;
    else
        y_pad = 0.08 * (y_max - y_min);
    end
    y_lims = [y_min - y_pad, y_max + y_pad];
end
