function out = plotTargetedDoGEndpointEffects(fig_dir, results, tests, contrast_lbl, precision_lbl, ps, varargin)
% plotTargetedDoGEndpointEffects  Dumbbell + forest summaries for targeted tests.
%
% Each metric gets one figure. The left panel shows absolute endpoint fitted
% values (L1 high information and L3 low information). The right panel shows
% the tested L3-L1 delta and CI from computeTargetedDoGHypothesisTests.

    ip = inputParser;
    addParameter(ip, 'n_back', NaN, @(x) isnumeric(x) && isscalar(x));
    addParameter(ip, 'show_middle_level', false, @(x) islogical(x) && isscalar(x));
    parse(ip, varargin{:});

    if nargin < 3 || isempty(tests)
        tests = computeTargetedDoGHypothesisTests(results, ip.Results.n_back);
    end
    if isempty(tests) || height(tests) == 0 || ~isfield(results, 'summary_table') || isempty(results.summary_table)
        out = struct('amplitude_pdf', '', 'fwhm_pdf', '');
        return
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    pack = packUnbinnedScatterParams(results.summary_table);
    plot_opts = local_buildPlotOpts(ps);
    ci_label = local_ciLabelFromTests(tests);

    out = struct();
    out.amplitude_pdf = local_renderMetric(fig_dir, pack.amp, tests, 'A', ...
        'Amplitude A (deg)', 'targeted_dog_endpoint_effects_amplitude.pdf', ...
        contrast_lbl, precision_lbl, plot_opts, ip.Results.n_back, ip.Results.show_middle_level, ci_label);
    out.fwhm_pdf = local_renderMetric(fig_dir, pack.fwhm, tests, 'FWHM', ...
        'FWHM (deg)', 'targeted_dog_endpoint_effects_fwhm.pdf', ...
        contrast_lbl, precision_lbl, plot_opts, ip.Results.n_back, ip.Results.show_middle_level, ci_label);
end

function out_pdf = local_renderMetric(fig_dir, vals, tests, parameter, x_label, fname, ...
    contrast_lbl, precision_lbl, plot_opts, n_back, show_middle_level, ci_label)

    axes_to_plot = {'previous_level', 'current_level', 'average_prev_curr'};
    axis_titles = {'Previous level endpoint', 'Current level endpoint', 'Average prev/current endpoint'};

    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 8.5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 8.5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 8.5]);
    tl = tiledlayout(numel(axes_to_plot), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    abs_vals = local_collectAbsoluteValues(vals, axes_to_plot, show_middle_level);
    delta_vals = local_collectDeltaValues(tests, parameter);
    abs_lims = local_limits(abs_vals, false);
    delta_lims = local_limits(delta_vals, true);

    for i_axis = 1:numel(axes_to_plot)
        axis_name = axes_to_plot{i_axis};
        ax_abs = nexttile(tl);
        local_plotDumbbell(ax_abs, vals, axis_name, abs_lims, x_label, ...
            contrast_lbl, precision_lbl, plot_opts, show_middle_level);
        title(ax_abs, axis_titles{i_axis}, 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');

        ax_delta = nexttile(tl);
        local_plotForest(ax_delta, tests, parameter, axis_name, delta_lims, plot_opts);
        title(ax_delta, 'Tested delta and CI', 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
    end

    if isfinite(n_back)
        fig_title = sprintf('%d-back targeted %s endpoint effects (%s)', n_back, parameter, ci_label);
    else
        fig_title = sprintf('Targeted %s endpoint effects (%s)', parameter, ci_label);
    end
    title(tl, fig_title, 'FontSize', plot_opts.axes_label_font_size + 1, 'Interpreter', 'none');

    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function local_plotDumbbell(ax, vals, axis_name, x_lims, x_label, contrast_lbl, precision_lbl, plot_opts, show_middle_level)
    hold(ax, 'on');
    y = [2, 1];
    labels = {'Contrast', 'Precision'};
    colors = [plot_opts.colors.blue; plot_opts.colors.green];
    level_labels = {contrast_lbl, precision_lbl};

    for m = 1:2
        v1 = local_axisLevelValue(vals(:, :, m), axis_name, 1);
        v2 = local_axisLevelValue(vals(:, :, m), axis_name, 2);
        v3 = local_axisLevelValue(vals(:, :, m), axis_name, 3);

        plot(ax, [v1, v3], [y(m), y(m)], '-', 'Color', colors(m, :) * 0.85, ...
            'LineWidth', plot_opts.line_width + 1, 'HandleVisibility', 'off');
        scatter(ax, v1, y(m), 60, 'MarkerFaceColor', colors(m, :), ...
            'MarkerEdgeColor', [1 1 1], 'LineWidth', plot_opts.line_width, ...
            'DisplayName', sprintf('%s high (%s)', labels{m}, level_labels{m}{1}));
        scatter(ax, v3, y(m), 60, 'MarkerFaceColor', colors(m, :) * 0.45, ...
            'MarkerEdgeColor', [1 1 1], 'LineWidth', plot_opts.line_width, ...
            'DisplayName', sprintf('%s low (%s)', labels{m}, level_labels{m}{3}));
        if show_middle_level
            scatter(ax, v2, y(m), 45, 'MarkerFaceColor', colors(m, :) * 0.65 + 0.35, ...
                'MarkerEdgeColor', colors(m, :) * 0.9, 'MarkerFaceAlpha', 0.35, ...
                'LineWidth', plot_opts.line_width, ...
                'DisplayName', sprintf('%s middle (%s)', labels{m}, level_labels{m}{2}));
        end
    end

    set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(labels));
    xlabel(ax, x_label, 'FontSize', plot_opts.axes_label_font_size);
    xlim(ax, x_lims);
    ylim(ax, [0.4, 2.6]);
    local_styleAxis(ax, plot_opts);
end

function local_plotForest(ax, tests, parameter, axis_name, x_lims, plot_opts)
    hold(ax, 'on');
    y = [5, 4, 3, 2, 1];
    row_names = local_effectNames(parameter);
    labels = local_effectLabels(parameter);
    colors = local_forestColors(plot_opts, numel(row_names));

    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    for i = 1:numel(row_names)
        row = local_getTestRow(tests, row_names{i}, parameter, axis_name);
        if isempty(row)
            continue
        end
        est = row.estimate(1);
        lo = row.bca_lo(1);
        hi = row.bca_hi(1);
        ok = local_isValid(row);
        if ok
            marker_face = colors(i, :);
            line_style = '-';
        else
            marker_face = [1 1 1];
            line_style = '--';
        end
        if isfinite(lo) && isfinite(hi)
            plot(ax, [lo, hi], [y(i), y(i)], line_style, ...
                'Color', colors(i, :) * 0.85, 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
        end
        if isfinite(est)
            scatter(ax, est, y(i), 45, 'MarkerFaceColor', marker_face, ...
                'MarkerEdgeColor', colors(i, :), 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
            sig_label = local_sigLabel(row);
            if ~isempty(sig_label)
                dx = 0.02 * diff(x_lims);
                text(ax, est + dx, y(i) + 0.18, sig_label, 'FontSize', 8, ...
                    'Color', colors(i, :) * 0.75, 'HorizontalAlignment', 'left', ...
                    'Interpreter', 'none');
            end
        end
        if ~ok
            if isfinite(hi)
                flag_x = hi;
            elseif isfinite(est)
                flag_x = est;
            else
                flag_x = 0;
            end
            text(ax, flag_x, y(i) + 0.18, 'flagged', 'FontSize', 8, ...
                'Color', plot_opts.colors.gray, 'HorizontalAlignment', 'right');
        end
    end

    set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(labels));
    xlabel(ax, sprintf('Delta %s (L3 - L1)', parameter), 'FontSize', plot_opts.axes_label_font_size, 'Interpreter', 'none');
    xlim(ax, x_lims);
    ylim(ax, [0.4, 5.6]);
    local_styleAxis(ax, plot_opts);
end

function v = local_axisLevelValue(mat, axis_name, level)
    switch axis_name
        case 'previous_level'
            v = mean(mat(level, :), 'omitnan');
        case 'current_level'
            v = mean(mat(:, level), 'omitnan');
        case 'average_prev_curr'
            v = 0.5 .* (mean(mat(level, :), 'omitnan') + mean(mat(:, level), 'omitnan'));
        otherwise
            v = NaN;
    end
end

function vals_out = local_collectAbsoluteValues(vals, axes_to_plot, show_middle_level)
    levels = [1, 3];
    if show_middle_level
        levels = [1, 2, 3];
    end
    vals_out = [];
    for i_axis = 1:numel(axes_to_plot)
        for m = 1:2
            for i_level = 1:numel(levels)
                vals_out(end+1, 1) = local_axisLevelValue(vals(:, :, m), axes_to_plot{i_axis}, levels(i_level)); %#ok<AGROW>
            end
        end
    end
end

function vals_out = local_collectDeltaValues(tests, parameter)
    vals_out = [];
    keep = strcmp(cellstr(tests.parameter), parameter);
    rows = tests(keep, :);
    for i = 1:height(rows)
        vals_out = [vals_out; rows.estimate(i); rows.bca_lo(i); rows.bca_hi(i)]; %#ok<AGROW>
    end
end

function row = local_getTestRow(tests, effect_name, parameter, axis_name)
    if contains(effect_name, 'boundary_interaction')
        target = sprintf('%s.%s', effect_name, 'boundary_interaction');
    else
        target = sprintf('%s.%s', effect_name, axis_name);
    end
    keep = strcmp(cellstr(tests.name), target) & strcmp(cellstr(tests.parameter), parameter);
    if any(keep)
        row = tests(find(keep, 1), :);
    else
        row = table();
    end
end

function names = local_effectNames(parameter)
    switch parameter
        case 'A'
            names = {'contrast_A_high_to_low', 'precision_A_high_to_low', ...
                'contrast_minus_precision_A', 'contrast_A_boundary_interaction', ...
                'precision_A_boundary_interaction'};
        case 'FWHM'
            names = {'contrast_FWHM_high_to_low', 'precision_FWHM_high_to_low', ...
                'precision_minus_contrast_FWHM', 'contrast_FWHM_boundary_interaction', ...
                'precision_FWHM_boundary_interaction'};
        otherwise
            names = {};
    end
end

function labels = local_effectLabels(parameter)
    switch parameter
        case 'A'
            labels = {'Contrast', 'Precision', 'Contrast - precision', ...
                'Contrast boundary', 'Precision boundary'};
        case 'FWHM'
            labels = {'Contrast', 'Precision', 'Precision - contrast', ...
                'Contrast boundary', 'Precision boundary'};
        otherwise
            labels = {};
    end
end

function colors = local_forestColors(plot_opts, n)
    base = [plot_opts.colors.blue; plot_opts.colors.green; plot_opts.colors.purple; ...
        plot_opts.colors.orange; plot_opts.colors.gray];
    if n <= size(base, 1)
        colors = base(1:n, :);
    else
        colors = repmat(base, ceil(n / size(base, 1)), 1);
        colors = colors(1:n, :);
    end
end

function tf = local_isValid(row)
    if ismember('valid_for_inference', row.Properties.VariableNames)
        tf = logical(row.valid_for_inference(1));
    else
        tf = true;
    end
end

function s = local_sigLabel(row)
    s = '';
    if ~local_isValid(row) || ~ismember('p_holm', row.Properties.VariableNames) || ...
            ~isfinite(row.p_holm(1)) || row.p_holm(1) >= 0.05
        return
    end
    if ismember('p_holm_label', row.Properties.VariableNames)
        s = sprintf('Holm p %s', char(row.p_holm_label(1)));
    else
        s = sprintf('Holm p %.3f', row.p_holm(1));
    end
end

function lims = local_limits(vals, include_zero)
    vals = vals(isfinite(vals));
    if isempty(vals)
        lims = [-1 1];
        return
    end
    if include_zero
        vals = [vals; 0];
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmin == vmax
        pad = max(abs(vmin), 1) * 0.1;
    else
        pad = 0.12 * (vmax - vmin);
    end
    lims = [vmin - pad, vmax + pad];
end

function local_styleAxis(ax, plot_opts)
    axis(ax, 'square');
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
        'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
end

function plot_opts = local_buildPlotOpts(ps)
    plot_opts.line_width = ps.line_width;
    plot_opts.colors = ps.colors;
    plot_opts.figure_color = ps.figure_color;
    plot_opts.tick_length = ps.tick_length;
    plot_opts.axes_label_font_size = 14;
    plot_opts.axes_tick_font_size = 13;
    if isfield(ps, 'axes_label_font_size')
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
    end
    if isfield(ps, 'axes_tick_font_size')
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    end
    if ~isfield(plot_opts.colors, 'purple')
        plot_opts.colors.purple = [102 51 204] ./ 255;
    end
    if ~isfield(plot_opts.colors, 'gray')
        plot_opts.colors.gray = [128 128 128] ./ 255;
    end
end

function ci_label = local_ciLabelFromTests(tests)
    ci_method = 'bootstrap';
    if istable(tests) && height(tests) > 0 && ismember('ci_method', tests.Properties.VariableNames)
        ci_method = tests.ci_method(1);
    end
    ci_label = sprintf('95%% %s CI', local_ciLabel(ci_method));
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
