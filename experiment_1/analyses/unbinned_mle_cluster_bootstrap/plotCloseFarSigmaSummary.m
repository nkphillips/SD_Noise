function plotCloseFarSigmaSummary(fig_dir, close_far_table, contrast_lbl, precision_lbl, ps, ci_pct, n_back, ci_method)
% plotCloseFarSigmaSummary  Close-vs-far sigma diagnostic per condition cell.

    if nargin < 6 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end
    if nargin < 7 || isempty(n_back)
        n_back = NaN;
    end
    if nargin < 8 || isempty(ci_method)
        ci_method = 'bootstrap';
    end
    ps = local_plotDefaults(ps);
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));
    fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 12 5], ...
        'PaperUnits', 'inches', 'PaperSize', [12 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 12 5]);

    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    manip_names = {'contrast', 'precision'};
    manip_titles = {'Contrast', 'Precision'};

    for m = 1:2
        ax = nexttile(tl);
        hold(ax, 'on');

        rows = strcmpi(cellstr(close_far_table.cond_manipulation), manip_names{m});
        Tm = close_far_table(rows, :);
        [~, ord] = sortrows([Tm.cond_prev, Tm.cond_curr], [1 2]);
        Tm = Tm(ord, :);

        x = 1:height(Tm);
        x_close = x - 0.12;
        x_far = x + 0.12;

        if m == 1
            base_color = ps.colors.blue;
            level_lbl = contrast_lbl;
        else
            base_color = ps.colors.green;
            level_lbl = precision_lbl;
        end
        close_color = 0.55 .* [1 1 1];
        far_color = base_color;

        for i = 1:height(Tm)
            plot(ax, [x_close(i), x_far(i)], [Tm.sigma_close(i), Tm.sigma_far(i)], '-', ...
                'Color', [0.75 0.75 0.75], 'LineWidth', ps.line_width, 'HandleVisibility', 'off');
        end

        local_errorScatter(ax, x_close, Tm.sigma_close, Tm.sigma_close_ci_lo, Tm.sigma_close_ci_hi, ...
            close_color, ps, 'Close');
        local_errorScatter(ax, x_far, Tm.sigma_far, Tm.sigma_far_ci_lo, Tm.sigma_far_ci_hi, ...
            far_color, ps, 'Far');

        yline(ax, 0, ':', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
        title(ax, manip_titles{m}, 'FontSize', ps.axes_tick_font_size);
        ylabel(ax, '\sigma (deg)', 'Interpreter', 'tex', 'FontSize', ps.axes_label_font_size);
        xlabel(ax, 'Previous-current level cell', 'FontSize', ps.axes_label_font_size);
        set(ax, 'XTick', x, 'XTickLabel', local_cellLabels(Tm, level_lbl), ...
            'XTickLabelRotation', 45);
        xlim(ax, [0.4, height(Tm) + 0.6]);
        local_setYLimits(ax, [Tm.sigma_close; Tm.sigma_far; Tm.sigma_close_ci_lo; ...
            Tm.sigma_close_ci_hi; Tm.sigma_far_ci_lo; Tm.sigma_far_ci_hi]);
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
            'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size);
        legend(ax, 'Location', 'best');
        text(ax, 0.02, 0.98, ci_label, ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Color', ps.colors.gray);
    end

    if isfinite(n_back)
        sgtitle(tl, sprintf('Close-vs-far psychometric \\sigma: %d-back (%s)', n_back, ci_label), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    else
        sgtitle(tl, sprintf('Close-vs-far psychometric \\sigma (%s)', ci_label), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    end

    exportgraphics(fig, fullfile(fig_dir, 'close_far_sigma_summary.pdf'), 'ContentType', 'vector');
    close(fig);

    local_plotDelta(fig_dir, close_far_table, contrast_lbl, precision_lbl, ps, ci_pct, n_back, ci_method);
    local_plotManipulationAverages(fig_dir, close_far_table, ps, n_back);
    local_plotManipulationDeltaAverages(fig_dir, close_far_table, ps, n_back);
    local_writeDeltaTestSummary(fig_dir, close_far_table);
end

function ps = local_plotDefaults(ps)
    if ~isfield(ps, 'figure_color') || isempty(ps.figure_color)
        ps.figure_color = [1 1 1];
    end
    if ~isfield(ps, 'line_width') || isempty(ps.line_width)
        ps.line_width = 1;
    end
    if ~isfield(ps, 'tick_length') || isempty(ps.tick_length)
        ps.tick_length = 0.020;
    end
    if ~isfield(ps, 'axes_label_font_size') || isempty(ps.axes_label_font_size)
        ps.axes_label_font_size = 14;
    end
    if ~isfield(ps, 'axes_tick_font_size') || isempty(ps.axes_tick_font_size)
        ps.axes_tick_font_size = 13;
    end
    if ~isfield(ps, 'colors') || isempty(ps.colors)
        ps.colors = struct();
    end
    if ~isfield(ps.colors, 'blue') || isempty(ps.colors.blue)
        ps.colors.blue = [38 71 237] / 255;
    end
    if ~isfield(ps.colors, 'green') || isempty(ps.colors.green)
        ps.colors.green = [0 153 0] / 255;
    end
    if ~isfield(ps.colors, 'gray') || isempty(ps.colors.gray)
        ps.colors.gray = [128 128 128] / 255;
    end
end

function local_errorScatter(ax, x, y, lo, hi, color, ps, label_str)
    yneg = max(0, y - lo);
    ypos = max(0, hi - y);
    errorbar(ax, x, y, yneg, ypos, 'o', ...
        'Color', color .* 0.8, 'MarkerFaceColor', color, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerSize', 5, 'LineWidth', ps.line_width, 'CapSize', 0, ...
        'DisplayName', label_str);
end

function labels = local_cellLabels(Tm, level_lbl)
    labels = cell(height(Tm), 1);
    for i = 1:height(Tm)
        labels{i} = sprintf('%s-%s', level_lbl{Tm.cond_prev(i)}, level_lbl{Tm.cond_curr(i)});
    end
end

function local_setYLimits(ax, values)
    values = values(isfinite(values));
    if isempty(values)
        ylim(ax, [0, 1]);
        return
    end
    y_min = min(values);
    y_max = max(values);
    if y_min == y_max
        pad = max(abs(y_min), 1) * 0.1;
    else
        pad = 0.08 * (y_max - y_min);
    end
    ylim(ax, [max(0, y_min - pad), y_max + pad]);
end

function local_plotDelta(fig_dir, close_far_table, contrast_lbl, precision_lbl, ps, ci_pct, n_back, ci_method)
    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));
    fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 12 5], ...
        'PaperUnits', 'inches', 'PaperSize', [12 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 12 5]);

    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    manip_names = {'contrast', 'precision'};
    manip_titles = {'Contrast', 'Precision'};
    for m = 1:2
        ax = nexttile(tl);
        hold(ax, 'on');

        rows = strcmpi(cellstr(close_far_table.cond_manipulation), manip_names{m});
        Tm = close_far_table(rows, :);
        [~, ord] = sortrows([Tm.cond_prev, Tm.cond_curr], [1 2]);
        Tm = Tm(ord, :);
        x = 1:height(Tm);

        if m == 1
            color = ps.colors.blue;
            level_lbl = contrast_lbl;
        else
            color = ps.colors.green;
            level_lbl = precision_lbl;
        end

        y = Tm.delta_sigma_far_minus_close;
        yneg = max(0, y - Tm.delta_sigma_ci_lo);
        ypos = max(0, Tm.delta_sigma_ci_hi - y);
        errorbar(ax, x, y, yneg, ypos, 'o', 'Color', color .* 0.8, ...
            'MarkerFaceColor', color, 'MarkerEdgeColor', [1 1 1], ...
            'MarkerSize', 5, 'LineWidth', ps.line_width, 'CapSize', 0);
        yline(ax, 0, ':', 'Color', [0.45 0.45 0.45]);

        title(ax, manip_titles{m}, 'FontSize', ps.axes_tick_font_size);
        ylabel(ax, '\Delta\sigma = \sigma_{far} - \sigma_{close} (deg)', ...
            'Interpreter', 'tex', 'FontSize', ps.axes_label_font_size);
        xlabel(ax, 'Previous-current level cell', 'FontSize', ps.axes_label_font_size);
        set(ax, 'XTick', x, 'XTickLabel', local_cellLabels(Tm, level_lbl), ...
            'XTickLabelRotation', 45);
        xlim(ax, [0.4, height(Tm) + 0.6]);
        local_setDeltaYLimits(ax, [Tm.delta_sigma_far_minus_close; Tm.delta_sigma_ci_lo; Tm.delta_sigma_ci_hi]);
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
            'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size);
        text(ax, 0.02, 0.98, ci_label, ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 7, 'Color', ps.colors.gray);
    end

    if isfinite(n_back)
        sgtitle(tl, sprintf('Close-far psychometric \\sigma difference: %d-back (%s)', n_back, ci_label), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    else
        sgtitle(tl, sprintf('Close-far psychometric \\sigma difference (%s)', ci_label), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    end

    exportgraphics(fig, fullfile(fig_dir, 'close_far_sigma_delta.pdf'), 'ContentType', 'vector');
    close(fig);
end

function local_setDeltaYLimits(ax, values)
    values = values(isfinite(values));
    if isempty(values)
        ylim(ax, [-1, 1]);
        return
    end
    lim = max(abs(values));
    if lim == 0
        lim = 1;
    end
    ylim(ax, [-1.1 * lim, 1.1 * lim]);
end

function local_plotManipulationAverages(fig_dir, close_far_table, ps, n_back)
    avg_table = local_manipulationAverageTable(close_far_table);
    if isempty(avg_table) || height(avg_table) == 0
        return
    end

    writetable(avg_table, fullfile(fig_dir, 'close_far_sigma_manipulation_average.csv'));

    fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 7 5], ...
        'PaperUnits', 'inches', 'PaperSize', [7 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 7 5]);
    ax = axes(fig);
    hold(ax, 'on');

    x = 1:height(avg_table);
    x_close = x - 0.14;
    x_far = x + 0.14;
    close_color = 0.55 .* [1 1 1];
    far_colors = nan(height(avg_table), 3);
    for i = 1:height(avg_table)
        if strcmpi(char(avg_table.cond_manipulation(i)), 'contrast')
            far_colors(i, :) = ps.colors.blue;
        else
            far_colors(i, :) = ps.colors.green;
        end
        plot(ax, [x_close(i), x_far(i)], ...
            [avg_table.sigma_close_mean(i), avg_table.sigma_far_mean(i)], '-', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', ps.line_width, ...
            'HandleVisibility', 'off');
    end

    errorbar(ax, x_close, avg_table.sigma_close_mean, avg_table.sigma_close_sem, ...
        'o', 'Color', close_color .* 0.8, 'MarkerFaceColor', close_color, ...
        'MarkerEdgeColor', [1 1 1], 'MarkerSize', 6, ...
        'LineWidth', ps.line_width, 'CapSize', 0, 'DisplayName', 'Close');
    for i = 1:height(avg_table)
        errorbar(ax, x_far(i), avg_table.sigma_far_mean(i), avg_table.sigma_far_sem(i), ...
            'o', 'Color', far_colors(i, :) .* 0.8, 'MarkerFaceColor', far_colors(i, :), ...
            'MarkerEdgeColor', [1 1 1], 'MarkerSize', 6, ...
            'LineWidth', ps.line_width, 'CapSize', 0, ...
            'HandleVisibility', local_handleVisibility(i));
    end

    y_vals = [avg_table.sigma_close_mean - avg_table.sigma_close_sem; ...
        avg_table.sigma_close_mean + avg_table.sigma_close_sem; ...
        avg_table.sigma_far_mean - avg_table.sigma_far_sem; ...
        avg_table.sigma_far_mean + avg_table.sigma_far_sem];
    local_setYLimits(ax, y_vals);
    y_lim = ylim(ax);
    y_text = y_lim(2) - 0.08 * diff(y_lim);
    for i = 1:height(avg_table)
        text(ax, x(i), y_text, sprintf('\\Delta\\sigma = %.2f \\pm %.2f', ...
            avg_table.delta_sigma_mean(i), avg_table.delta_sigma_sem(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Color', ps.colors.gray, 'Interpreter', 'tex');
    end

    set(ax, 'XTick', x, 'XTickLabel', local_titleCase(cellstr(string(avg_table.cond_manipulation))));
    ylabel(ax, 'Mean \sigma across cells (deg)', 'Interpreter', 'tex', ...
        'FontSize', ps.axes_label_font_size);
    xlabel(ax, 'Manipulation', 'FontSize', ps.axes_label_font_size);
    xlim(ax, [0.5, height(avg_table) + 0.5]);
    axis(ax, 'square');
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
        'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size);
    legend(ax, {'Close', 'Far'}, 'Location', 'best');

    if isfinite(n_back)
        title(ax, sprintf('Close-vs-far \\sigma averaged across cells: %d-back', n_back), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    else
        title(ax, 'Close-vs-far \sigma averaged across cells', ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    end

    exportgraphics(fig, fullfile(fig_dir, 'close_far_sigma_manipulation_average.pdf'), ...
        'ContentType', 'vector');
    close(fig);
end

function avg_table = local_manipulationAverageTable(close_far_table)
    manip_names = {'contrast', 'precision'};
    manip_col = cell(0, 1);
    n_cells = zeros(0, 1);
    sigma_close_mean = zeros(0, 1);
    sigma_close_sem = zeros(0, 1);
    sigma_far_mean = zeros(0, 1);
    sigma_far_sem = zeros(0, 1);
    delta_sigma_mean = zeros(0, 1);
    delta_sigma_sem = zeros(0, 1);
    delta_sigma_t = zeros(0, 1);
    delta_sigma_df = zeros(0, 1);
    delta_sigma_p = zeros(0, 1);
    delta_sigma_ci_lo = zeros(0, 1);
    delta_sigma_ci_hi = zeros(0, 1);

    for m = 1:numel(manip_names)
        rows = strcmpi(cellstr(close_far_table.cond_manipulation), manip_names{m});
        Tm = close_far_table(rows, :);
        if isempty(Tm) || height(Tm) == 0
            continue
        end

        [close_mean, close_sem] = local_meanSem(Tm.sigma_close);
        [far_mean, far_sem] = local_meanSem(Tm.sigma_far);
        [delta_mean, delta_sem, delta_n] = local_meanSem(Tm.delta_sigma_far_minus_close);
        delta_test = local_oneSampleTTest(Tm.delta_sigma_far_minus_close);

        manip_col{end+1, 1} = manip_names{m}; %#ok<AGROW>
        n_cells(end+1, 1) = delta_n; %#ok<AGROW>
        sigma_close_mean(end+1, 1) = close_mean; %#ok<AGROW>
        sigma_close_sem(end+1, 1) = close_sem; %#ok<AGROW>
        sigma_far_mean(end+1, 1) = far_mean; %#ok<AGROW>
        sigma_far_sem(end+1, 1) = far_sem; %#ok<AGROW>
        delta_sigma_mean(end+1, 1) = delta_mean; %#ok<AGROW>
        delta_sigma_sem(end+1, 1) = delta_sem; %#ok<AGROW>
        delta_sigma_t(end+1, 1) = delta_test.t; %#ok<AGROW>
        delta_sigma_df(end+1, 1) = delta_test.df; %#ok<AGROW>
        delta_sigma_p(end+1, 1) = delta_test.p; %#ok<AGROW>
        delta_sigma_ci_lo(end+1, 1) = delta_test.ci_lo; %#ok<AGROW>
        delta_sigma_ci_hi(end+1, 1) = delta_test.ci_hi; %#ok<AGROW>
    end

    cond_manipulation = categorical(manip_col);
    avg_table = table(cond_manipulation, n_cells, sigma_close_mean, sigma_close_sem, ...
        sigma_far_mean, sigma_far_sem, delta_sigma_mean, delta_sigma_sem, ...
        delta_sigma_t, delta_sigma_df, delta_sigma_p, delta_sigma_ci_lo, delta_sigma_ci_hi, ...
        'VariableNames', {'cond_manipulation', 'n_cells', ...
        'sigma_close_mean', 'sigma_close_sem', ...
        'sigma_far_mean', 'sigma_far_sem', ...
        'delta_sigma_mean', 'delta_sigma_sem', ...
        'delta_sigma_t', 'delta_sigma_df', 'delta_sigma_p', ...
        'delta_sigma_ci_lo', 'delta_sigma_ci_hi'});
end

function [mu, sem, n] = local_meanSem(values)
    values = values(isfinite(values));
    values = values(:);
    n = numel(values);
    if n == 0
        mu = NaN;
        sem = NaN;
    elseif n == 1
        mu = values;
        sem = NaN;
    else
        mu = mean(values);
        sem = std(values, 0) ./ sqrt(n);
    end
end

function stats = local_oneSampleTTest(values)
    values = values(isfinite(values));
    values = values(:);
    n = numel(values);
    stats = struct('t', NaN, 'df', max(n - 1, 0), 'p', NaN, ...
        'ci_lo', NaN, 'ci_hi', NaN);
    if n < 2
        return
    end

    mu = mean(values);
    sd = std(values, 0);
    sem = sd / sqrt(n);
    stats.df = n - 1;
    if sem == 0
        if mu > 0
            stats.t = Inf;
            stats.p = 0;
        elseif mu < 0
            stats.t = -Inf;
            stats.p = 0;
        else
            stats.t = 0;
            stats.p = 1;
        end
        stats.ci_lo = mu;
        stats.ci_hi = mu;
        return
    end

    stats.t = mu / sem;
    stats.p = 2 * (1 - tcdf(abs(stats.t), stats.df));
    tcrit = tinv(0.975, stats.df);
    stats.ci_lo = mu - tcrit * sem;
    stats.ci_hi = mu + tcrit * sem;
end

function out = local_handleVisibility(i)
    if i == 1
        out = 'on';
    else
        out = 'off';
    end
end

function labels = local_titleCase(labels)
    for i = 1:numel(labels)
        if isempty(labels{i})
            continue
        end
        labels{i} = [upper(labels{i}(1)), lower(labels{i}(2:end))];
    end
end

function local_plotManipulationDeltaAverages(fig_dir, close_far_table, ps, n_back)
    avg_table = local_manipulationAverageTable(close_far_table);
    if isempty(avg_table) || height(avg_table) == 0
        return
    end

    fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 9 4.8], ...
        'PaperUnits', 'inches', 'PaperSize', [9 4.8], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 9 4.8]);

    tl = tiledlayout(fig, 1, height(avg_table), 'TileSpacing', 'compact', 'Padding', 'compact');
    y_vals = [avg_table.delta_sigma_mean - avg_table.delta_sigma_sem; ...
        avg_table.delta_sigma_mean + avg_table.delta_sigma_sem; ...
        avg_table.delta_sigma_ci_lo; avg_table.delta_sigma_ci_hi; 0];
    y_vals = y_vals(isfinite(y_vals));
    if isempty(y_vals)
        y_lim = [-1, 1];
    else
        lim = max(abs(y_vals));
        if lim == 0
            lim = 1;
        end
        y_lim = [-1.35 * lim, 1.35 * lim];
    end

    for i = 1:height(avg_table)
        ax = nexttile(tl);
        hold(ax, 'on');

        manip_name = char(avg_table.cond_manipulation(i));
        if strcmpi(manip_name, 'contrast')
            color = ps.colors.blue;
        else
            color = ps.colors.green;
        end
        title_label = local_titleCase({manip_name});
        title_label = title_label{1};

        errorbar(ax, 1, avg_table.delta_sigma_mean(i), avg_table.delta_sigma_sem(i), ...
            'o', 'Color', color .* 0.8, 'MarkerFaceColor', color, ...
            'MarkerEdgeColor', [1 1 1], 'MarkerSize', 8, ...
            'LineWidth', ps.line_width, 'CapSize', 0);
        yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');

        ylim(ax, y_lim);
        xlim(ax, [0.55, 1.45]);
        set(ax, 'XTick', 1, 'XTickLabel', {title_label});
        ylabel(ax, '\Delta\sigma = \sigma_{far} - \sigma_{close} (deg)', ...
            'Interpreter', 'tex', 'FontSize', ps.axes_label_font_size);
        title(ax, title_label, 'FontSize', ps.axes_tick_font_size);
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
            'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size);

        p_str = local_formatP(avg_table.delta_sigma_p(i));
        stats_str = sprintf('mean ± SEM = %.2f ± %.2f\nt(%d) = %.2f, p = %s', ...
            avg_table.delta_sigma_mean(i), avg_table.delta_sigma_sem(i), ...
            round(avg_table.delta_sigma_df(i)), avg_table.delta_sigma_t(i), p_str);
        text(ax, 0.05, 0.95, stats_str, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'Color', ps.colors.gray, 'Interpreter', 'none', ...
            'BackgroundColor', ps.figure_color, 'Margin', 3);
    end

    if isfinite(n_back)
        sgtitle(tl, sprintf('Mean close-far \\sigma difference across cells: %d-back', n_back), ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    else
        sgtitle(tl, 'Mean close-far \sigma difference across cells', ...
            'FontSize', ps.axes_label_font_size + 1, 'Interpreter', 'tex');
    end

    exportgraphics(fig, fullfile(fig_dir, 'close_far_sigma_delta_manipulation_average.pdf'), ...
        'ContentType', 'vector');
    close(fig);
end

function p_str = local_formatP(p)
    if ~isfinite(p)
        p_str = 'NA';
    elseif p < 0.001
        p_str = '<0.001';
    else
        p_str = sprintf('%.3f', p);
    end
end

function local_writeDeltaTestSummary(fig_dir, close_far_table)
    test_table = local_deltaTestSummaryTable(close_far_table);
    if isempty(test_table) || height(test_table) == 0
        return
    end

    writetable(test_table, fullfile(fig_dir, 'close_far_sigma_delta_tests.csv'));
end

function test_table = local_deltaTestSummaryTable(close_far_table)
    manip_names = {'contrast', 'precision'};

    test_name = strings(0, 1);
    test_family = strings(0, 1);
    comparison = strings(0, 1);
    n_cells = zeros(0, 1);
    mean_delta_sigma = nan(0, 1);
    sem_delta_sigma = nan(0, 1);
    mean_difference = nan(0, 1);
    sem_difference = nan(0, 1);
    t_stat = nan(0, 1);
    df = nan(0, 1);
    p_value = nan(0, 1);
    ci_lo = nan(0, 1);
    ci_hi = nan(0, 1);

    for m = 1:numel(manip_names)
        vals = local_deltaSigmaForManip(close_far_table, manip_names{m});
        [mu, sem, n] = local_meanSem(vals);
        stats = local_oneSampleTTest(vals);

        test_name(end+1, 1) = sprintf('%s_delta_sigma_vs_zero', manip_names{m}); %#ok<AGROW>
        test_family(end+1, 1) = "within_manipulation"; %#ok<AGROW>
        comparison(end+1, 1) = sprintf('%s: delta_sigma = 0', manip_names{m}); %#ok<AGROW>
        n_cells(end+1, 1) = n; %#ok<AGROW>
        mean_delta_sigma(end+1, 1) = mu; %#ok<AGROW>
        sem_delta_sigma(end+1, 1) = sem; %#ok<AGROW>
        mean_difference(end+1, 1) = mu; %#ok<AGROW>
        sem_difference(end+1, 1) = sem; %#ok<AGROW>
        t_stat(end+1, 1) = stats.t; %#ok<AGROW>
        df(end+1, 1) = stats.df; %#ok<AGROW>
        p_value(end+1, 1) = stats.p; %#ok<AGROW>
        ci_lo(end+1, 1) = stats.ci_lo; %#ok<AGROW>
        ci_hi(end+1, 1) = stats.ci_hi; %#ok<AGROW>
    end

    [contrast_vals, precision_vals] = local_matchedDeltaSigma(close_far_table);
    diff_vals = precision_vals - contrast_vals;
    [mu_diff, sem_diff, n_diff] = local_meanSem(diff_vals);
    stats_diff = local_oneSampleTTest(diff_vals);

    test_name(end+1, 1) = "precision_minus_contrast_delta_sigma"; %#ok<AGROW>
    test_family(end+1, 1) = "between_manipulation"; %#ok<AGROW>
    comparison(end+1, 1) = "matched cell ranks: precision delta_sigma - contrast delta_sigma = 0"; %#ok<AGROW>
    n_cells(end+1, 1) = n_diff; %#ok<AGROW>
    mean_delta_sigma(end+1, 1) = NaN; %#ok<AGROW>
    sem_delta_sigma(end+1, 1) = NaN; %#ok<AGROW>
    mean_difference(end+1, 1) = mu_diff; %#ok<AGROW>
    sem_difference(end+1, 1) = sem_diff; %#ok<AGROW>
    t_stat(end+1, 1) = stats_diff.t; %#ok<AGROW>
    df(end+1, 1) = stats_diff.df; %#ok<AGROW>
    p_value(end+1, 1) = stats_diff.p; %#ok<AGROW>
    ci_lo(end+1, 1) = stats_diff.ci_lo; %#ok<AGROW>
    ci_hi(end+1, 1) = stats_diff.ci_hi; %#ok<AGROW>

    test_table = table(test_name, test_family, comparison, n_cells, ...
        mean_delta_sigma, sem_delta_sigma, mean_difference, sem_difference, ...
        t_stat, df, p_value, ci_lo, ci_hi);
end

function vals = local_deltaSigmaForManip(close_far_table, manip_name)
    rows = strcmpi(cellstr(close_far_table.cond_manipulation), manip_name);
    Tm = close_far_table(rows, :);
    if isempty(Tm) || height(Tm) == 0
        vals = [];
        return
    end

    [~, ord] = sortrows([Tm.cond_prev, Tm.cond_curr], [1 2]);
    Tm = Tm(ord, :);
    vals = Tm.delta_sigma_far_minus_close(:);
    vals = vals(isfinite(vals));
end

function [contrast_vals, precision_vals] = local_matchedDeltaSigma(close_far_table)
    Tc = local_sortedManipTable(close_far_table, 'contrast');
    Tp = local_sortedManipTable(close_far_table, 'precision');
    contrast_vals = [];
    precision_vals = [];
    if isempty(Tc) || isempty(Tp)
        return
    end

    for i = 1:height(Tc)
        match = Tp.cond_prev == Tc.cond_prev(i) & Tp.cond_curr == Tc.cond_curr(i);
        if ~any(match)
            continue
        end

        cv = Tc.delta_sigma_far_minus_close(i);
        pv = Tp.delta_sigma_far_minus_close(find(match, 1, 'first'));
        if isfinite(cv) && isfinite(pv)
            contrast_vals(end+1, 1) = cv; %#ok<AGROW>
            precision_vals(end+1, 1) = pv; %#ok<AGROW>
        end
    end
end

function Tm = local_sortedManipTable(close_far_table, manip_name)
    rows = strcmpi(cellstr(close_far_table.cond_manipulation), manip_name);
    Tm = close_far_table(rows, :);
    if isempty(Tm) || height(Tm) == 0
        return
    end
    [~, ord] = sortrows([Tm.cond_prev, Tm.cond_curr], [1 2]);
    Tm = Tm(ord, :);
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
