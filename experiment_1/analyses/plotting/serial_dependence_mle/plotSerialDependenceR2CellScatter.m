function plotSerialDependenceR2CellScatter(fig_dir, summary_table, ps, n_back)
% plotSerialDependenceR2CellScatter  Per-cell DoG R2 diagnostic scatter.
%
% Renders two panels, one per manipulation, with the 3 x 3 previous-current
% cells shown as nine points. Uses r2_delta_bins, the Efron R2 aligned with
% the plotted DoG-vs-delta-theta curves.

    if nargin < 4 || isempty(n_back)
        n_back = NaN;
    end
    if nargin < 3 || isempty(ps)
        ps = plotSettings();
    end
    if isempty(summary_table) || ~ismember('r2_delta_bins', summary_table.Properties.VariableNames)
        return
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    axes_fs = 14;
    tick_fs = 13;
    try
        axes_fs = ps.axes_label_font_size;
        tick_fs = ps.axes_tick_font_size;
    catch
    end

    fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 9 4.5], ...
        'PaperUnits', 'inches', 'PaperSize', [9 4.5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 9 4.5]);
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    y_all = summary_table.r2_delta_bins;
    y_all = y_all(isfinite(y_all));
    if isempty(y_all)
        y_lims = [0, 1];
    else
        y_min = min(y_all);
        y_pad = 0.05 * max(1 - y_min, 0.01);
        y_lims = [max(0, y_min - y_pad), 1];
    end

    manip_names = {'contrast', 'precision'};
    manip_titles = {'Contrast', 'Precision'};
    manip_colors = {ps.colors.blue, ps.colors.green};
    x_labels = local_cellLabels();

    for m = 1:2
        ax = nexttile(tl);
        hold(ax, 'on');

        keep = strcmpi(cellstr(string(summary_table.cond_manipulation)), manip_names{m});
        sub = summary_table(keep, :);
        [x, y] = local_orderedR2(sub);
        ok = isfinite(y);

        scatter(ax, x(ok), y(ok), 55, ...
            'MarkerFaceColor', manip_colors{m}, ...
            'MarkerEdgeColor', [1 1 1], ...
            'MarkerFaceAlpha', 0.80, ...
            'LineWidth', ps.line_width);

        vals = y(ok);
        if ~isempty(vals)
            mean_r2 = mean(vals);
            sd_r2 = std(vals, 0);
            sem_r2 = sd_r2 / sqrt(numel(vals));
            yline(ax, mean_r2, '-', 'Color', ps.colors.black, ...
                'LineWidth', ps.line_width, 'HandleVisibility', 'off');
            yline(ax, mean_r2 - sem_r2, ':', 'Color', ps.colors.gray, ...
                'LineWidth', ps.line_width, 'HandleVisibility', 'off');
            yline(ax, mean_r2 + sem_r2, ':', 'Color', ps.colors.gray, ...
                'LineWidth', ps.line_width, 'HandleVisibility', 'off');
            [txt_x, txt_y, txt_va] = local_annotationPosition(mean_r2, sem_r2, y_lims);
            text(ax, txt_x, txt_y, sprintf('mean = %.2f%sSEM = %.2f', mean_r2, newline, sem_r2), ...
                'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', txt_va, 'FontSize', 9, ...
                'Color', ps.colors.black, 'Interpreter', 'none', ...
                'BackgroundColor', ps.figure_color, 'Margin', 2);
        end

        xlim(ax, [0.5, 9.5]);
        ylim(ax, y_lims);
        set(ax, 'XTick', 1:9, 'XTickLabel', x_labels, ...
            'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
            'LineWidth', ps.line_width, 'FontSize', tick_fs);
        xtickangle(ax, 45);
        xlabel(ax, 'Previous-current cell', 'FontSize', axes_fs);
        if m == 1
            ylabel(ax, 'DoG R^2_{\Delta\theta}', 'FontSize', axes_fs, 'Interpreter', 'tex');
        end
        title(ax, manip_titles{m}, 'FontSize', tick_fs);
        axis(ax, 'square');
        box(ax, 'off');
    end

    if isfinite(n_back)
        title_str = sprintf('Serial dependence fit quality: %d-back', n_back);
    else
        title_str = 'Serial dependence fit quality';
    end
    title(tl, title_str, 'FontSize', axes_fs + 1);

    out_pdf = fullfile(fig_dir, 'serial_dependence_fit_quality_cell_scatter.pdf');
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function labels = local_cellLabels()
    labels = cell(1, 9);
    k = 0;
    for prev = 1:3
        for curr = 1:3
            k = k + 1;
            labels{k} = sprintf('P%d-C%d', prev, curr);
        end
    end
end

function [x, y] = local_orderedR2(sub)
    x = 1:9;
    y = nan(1, 9);
    k = 0;
    for prev = 1:3
        for curr = 1:3
            k = k + 1;
            idx = sub.cond_prev == prev & sub.cond_curr == curr;
            if any(idx)
                vals = sub.r2_delta_bins(idx);
                y(k) = vals(1);
            end
        end
    end
end

function [txt_x, txt_y, txt_va] = local_annotationPosition(mean_r2, sem_r2, y_lims)
    txt_x = 0.04;
    band_lo = mean_r2 - sem_r2;
    band_hi = mean_r2 + sem_r2;
    y_range = max(diff(y_lims), eps);
    top_clearance = (y_lims(2) - band_hi) / y_range;
    bottom_clearance = (band_lo - y_lims(1)) / y_range;

    if top_clearance >= bottom_clearance
        txt_y = 0.96;
        txt_va = 'top';
    else
        txt_y = 0.04;
        txt_va = 'bottom';
    end
end
