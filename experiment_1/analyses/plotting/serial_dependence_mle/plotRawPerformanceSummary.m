function out = plotRawPerformanceSummary(perf, fig_dir, ps)
% plotRawPerformanceSummary  Save raw performance CSVs and summary figures.

    if nargin < 3 || isempty(ps)
        ps = plotSettings();
    end
    ps = local_completePlotSettings(ps);
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    out = struct();
    out.csvs = local_writeCsvs(perf, fig_dir);
    out.overall_pdf = local_plotOverall(perf, fig_dir, ps);
    out.comparison_pdf = local_plotManipulationComparison(perf, fig_dir, ps);
    out.group_level_pair_pdf = local_plotLevelPairHeatmaps(perf.group_level_pair, fig_dir, ps, ...
        'raw_performance_level_pair_group.pdf', 'Group mean percent correct', 'mean_accuracy', 'total_accuracy_trials');
    out.super_subject_level_pair_pdf = local_plotLevelPairHeatmaps(perf.super_subject_level_pair, fig_dir, ps, ...
        'raw_performance_level_pair_super_subject.pdf', 'Super-subject percent correct', 'accuracy', 'n_accuracy_trials');
end

function csvs = local_writeCsvs(perf, fig_dir)
    table_fields = fieldnames(perf);
    csvs = struct();
    for i = 1:numel(table_fields)
        field = table_fields{i};
        if istable(perf.(field))
            csv_path = fullfile(fig_dir, sprintf('raw_performance_%s.csv', field));
            writetable(perf.(field), csv_path);
            csvs.(field) = csv_path;
        end
    end
end

function ps = local_completePlotSettings(ps)
    if ~isfield(ps, 'font_type') || isempty(ps.font_type)
        ps.font_type = 'Helvetica';
    end
    if ~isfield(ps, 'axes_label_font_size') || isempty(ps.axes_label_font_size)
        ps.axes_label_font_size = 14;
    end
    if ~isfield(ps, 'axes_tick_font_size') || isempty(ps.axes_tick_font_size)
        ps.axes_tick_font_size = 13;
    end
    if ~isfield(ps, 'tick_length') || isempty(ps.tick_length)
        ps.tick_length = 0.020;
    end
    if ~isfield(ps, 'line_width') || isempty(ps.line_width)
        ps.line_width = 1;
    end
    if ~isfield(ps, 'colors') || ~isstruct(ps.colors)
        ps.colors = struct();
    end
    ps.colors = local_completeColors(ps.colors);
end

function colors = local_completeColors(colors)
    colors = local_defaultColor(colors, 'blue', [0 76 152] / 255);
    colors = local_defaultColor(colors, 'green', [0 153 0] / 255);
    colors = local_defaultColor(colors, 'white', [1 1 1]);
end

function colors = local_defaultColor(colors, field, value)
    if ~isfield(colors, field) || isempty(colors.(field))
        colors.(field) = value;
    end
end

function pdf_path = local_plotOverall(perf, fig_dir, ps)
    subj_overall = perf.subject_overall;
    subj_manip = perf.subject_manipulation;
    group_overall = perf.group_overall;
    group_manip = perf.group_manipulation;
    subj_ids = subj_overall.subject_id;
    x_labels = {'Overall', 'Contrast', 'Precision'};
    colors = [0.35 0.35 0.35; ps.colors.blue; ps.colors.green];

    y = nan(numel(subj_ids), 3);
    y(:, 1) = subj_overall.accuracy;
    for s = 1:numel(subj_ids)
        c_row = subj_manip(subj_manip.subject_id == subj_ids(s) & subj_manip.cond_manipulation == 'contrast', :);
        p_row = subj_manip(subj_manip.subject_id == subj_ids(s) & subj_manip.cond_manipulation == 'precision', :);
        if ~isempty(c_row); y(s, 2) = c_row.accuracy(1); end
        if ~isempty(p_row); y(s, 3) = p_row.accuracy(1); end
    end

    group_mean = nan(1, 3);
    group_sem = nan(1, 3);
    group_mean(1) = group_overall.mean_accuracy(1);
    group_sem(1) = group_overall.sem_accuracy(1);
    for m = 1:2
        manip = {'contrast', 'precision'};
        row = group_manip(group_manip.cond_manipulation == manip{m}, :);
        if ~isempty(row)
            group_mean(m + 1) = row.mean_accuracy(1);
            group_sem(m + 1) = row.sem_accuracy(1);
        end
    end

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 7.2 5.8], ...
        'PaperUnits', 'inches', 'PaperSize', [7.2 5.8], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 7.2 5.8]);
    ax = axes(fig);
    hold(ax, 'on');
    for x = 1:3
        jitter = linspace(-0.10, 0.10, max(numel(subj_ids), 2))';
        jitter = jitter(1:numel(subj_ids));
        scatter(ax, x + jitter, 100 * y(:, x), 28, ...
            'MarkerFaceColor', [0.72 0.72 0.72], ...
            'MarkerEdgeColor', [0.30 0.30 0.30], ...
            'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.65);
        errorbar(ax, x, 100 * group_mean(x), 100 * group_sem(x), 'o', ...
            'Color', colors(x, :), ...
            'MarkerFaceColor', colors(x, :), ...
            'MarkerEdgeColor', ps.colors.white, ...
            'LineWidth', 1.4, 'MarkerSize', 8, 'CapSize', 0);
    end
    yline(ax, 50, 'k:', 'HandleVisibility', 'off');
    xlim(ax, [0.5 3.5]);
    ylim(ax, [0 100]);
    set(ax, 'XTick', 1:3, 'XTickLabel', x_labels);
    ylabel(ax, 'Percent correct', 'FontSize', ps.axes_label_font_size);
    title(ax, 'Raw 2AFC performance', 'FontSize', ps.axes_label_font_size);
    local_styleAxis(ax, ps);

    pdf_path = fullfile(fig_dir, 'raw_performance_overall.pdf');
    local_exportFigure(fig, pdf_path);
end

function pdf_path = local_plotManipulationComparison(perf, fig_dir, ps)
    subj_manip = perf.subject_manipulation;
    comp = perf.manipulation_comparison;
    subj_ids = unique(subj_manip.subject_id, 'stable');
    y_contrast = nan(numel(subj_ids), 1);
    y_precision = nan(numel(subj_ids), 1);
    for s = 1:numel(subj_ids)
        c_row = subj_manip(subj_manip.subject_id == subj_ids(s) & subj_manip.cond_manipulation == 'contrast', :);
        p_row = subj_manip(subj_manip.subject_id == subj_ids(s) & subj_manip.cond_manipulation == 'precision', :);
        if ~isempty(c_row); y_contrast(s) = c_row.accuracy(1); end
        if ~isempty(p_row); y_precision(s) = p_row.accuracy(1); end
    end
    valid = isfinite(y_contrast) & isfinite(y_precision);

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 6.2 5.8], ...
        'PaperUnits', 'inches', 'PaperSize', [6.2 5.8], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 6.2 5.8]);
    ax = axes(fig);
    hold(ax, 'on');
    for s = find(valid)'
        plot(ax, [1 2], 100 * [y_contrast(s), y_precision(s)], '-', ...
            'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
    end
    scatter(ax, ones(sum(valid), 1), 100 * y_contrast(valid), 34, ...
        'MarkerFaceColor', ps.colors.blue, 'MarkerEdgeColor', ps.colors.white);
    scatter(ax, 2 * ones(sum(valid), 1), 100 * y_precision(valid), 34, ...
        'MarkerFaceColor', ps.colors.green, 'MarkerEdgeColor', ps.colors.white);

    if ~isempty(comp) && isfinite(comp.mean_difference(1))
        txt = sprintf('Precision - contrast: %.1f pp, 95%% CI [%.1f, %.1f], p = %s', ...
            100 * comp.mean_difference(1), 100 * comp.ci95_lo(1), 100 * comp.ci95_hi(1), ...
            local_formatP(comp.p_value(1)));
        text(ax, 0.5, 0.96, txt, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', ps.axes_tick_font_size - 1, 'Interpreter', 'none');
    end

    yline(ax, 50, 'k:', 'HandleVisibility', 'off');
    xlim(ax, [0.6 2.4]);
    ylim(ax, [0 100]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Contrast', 'Precision'});
    ylabel(ax, 'Percent correct', 'FontSize', ps.axes_label_font_size);
    title(ax, 'Matched-subject performance comparison', 'FontSize', ps.axes_label_font_size);
    local_styleAxis(ax, ps);

    pdf_path = fullfile(fig_dir, 'raw_performance_manipulation_comparison.pdf');
    local_exportFigure(fig, pdf_path);
end

function pdf_path = local_plotLevelPairHeatmaps(T, fig_dir, ps, file_name, title_prefix, value_var, n_var)
    manip_names = {'contrast', 'precision'};
    manip_labels = {'Contrast', 'Precision'};

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 10.5 4.8], ...
        'PaperUnits', 'inches', 'PaperSize', [10.5 4.8], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 10.5 4.8]);
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    for m = 1:2
        ax = nexttile(tl);
        mat = nan(3, 3);
        n_mat = nan(3, 3);
        for prev_lvl = 1:3
            for curr_lvl = 1:3
                row = T(T.cond_manipulation == manip_names{m} & ...
                    T.cond_prev == prev_lvl & T.cond_curr == curr_lvl, :);
                if ~isempty(row)
                    mat(prev_lvl, curr_lvl) = row.(value_var)(1);
                    n_mat(prev_lvl, curr_lvl) = row.(n_var)(1);
                end
            end
        end

        imagesc(ax, 100 * mat);
        colormap(ax, local_blueWhiteCmap(ps.colors.blue));
        caxis(ax, [50 100]);
        axis(ax, 'square');
        set(ax, 'XTick', 1:3, 'XTickLabel', {'L1', 'L2', 'L3'}, ...
            'YTick', 1:3, 'YTickLabel', {'L1', 'L2', 'L3'});
        xlabel(ax, 'Current level', 'FontSize', ps.axes_label_font_size);
        ylabel(ax, 'Previous level', 'FontSize', ps.axes_label_font_size);
        title(ax, sprintf('%s: %s', title_prefix, manip_labels{m}), ...
            'FontSize', ps.axes_label_font_size, 'Interpreter', 'none');
        for prev_lvl = 1:3
            for curr_lvl = 1:3
                if isfinite(mat(prev_lvl, curr_lvl))
                    text(ax, curr_lvl, prev_lvl, sprintf('%.1f%%\nn=%d', ...
                        100 * mat(prev_lvl, curr_lvl), round(n_mat(prev_lvl, curr_lvl))), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', ps.axes_tick_font_size - 2, ...
                        'Color', local_textColor(mat(prev_lvl, curr_lvl)));
                end
            end
        end
        cb = colorbar(ax);
        cb.Label.String = 'Percent correct';
        cb.Label.FontSize = ps.axes_label_font_size;
        local_styleAxis(ax, ps);
    end

    title(tl, sprintf('%s by level pair', title_prefix), ...
        'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'none');
    pdf_path = fullfile(fig_dir, file_name);
    local_exportFigure(fig, pdf_path);
end

function cmap = local_blueWhiteCmap(blue)
    n = 256;
    white = [1 1 1];
    cmap = [linspace(white(1), blue(1), n)', ...
            linspace(white(2), blue(2), n)', ...
            linspace(white(3), blue(3), n)'];
end

function color = local_textColor(acc)
    if acc >= 0.75
        color = [1 1 1];
    else
        color = [0 0 0];
    end
end

function p = local_formatP(p_value)
    if ~isfinite(p_value)
        p = 'NA';
    elseif p_value < 0.001
        p = '< .001';
    else
        p = sprintf('%.3f', p_value);
    end
end

function local_styleAxis(ax, ps)
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
        'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size, ...
        'FontName', ps.font_type);
end

function local_exportFigure(fig, pdf_path)
    try
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    catch exportErr
        warning('plotRawPerformanceSummary:vectorExportFailed', ...
            'Vector export failed (%s). Falling back to raster PDF.', exportErr.message);
        exportgraphics(fig, pdf_path, 'ContentType', 'image', 'Resolution', 200);
    end
    close(fig);
end
