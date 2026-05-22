function plotSubjectDataQuality(tbl_trials, fig_dir, subj_labels, baselines)
% plotSubjectDataQuality  Per-subject data quality diagnostic.
%
% For each subject, generates one row with three panels:
%   (1) Trial counts per (manipulation, prev, curr) cell -- 18 bars, blue=contrast,
%       green=precision. Identifies sparse cells.
%   (2) Empirical psychometric P(CW) vs binned x_probe, with marker size scaled
%       by trial count at that x. If a baseline fit is supplied, the fitted
%       Gaussian psychometric (with 25% lapse) is overlaid in red and the
%       baseline mu_i is marked.
%   (3) Summary statistics: total trials, response rate, accuracy, median /
%       minimum trials per cell, baseline mu_i and sigma_i, skipped flag.
%
% Saves a single multi-row PDF and a per-subject CSV summary. Designed to be
% read top-down so you can spot anomalies (sparse data, near-chance accuracy,
% flat or non-monotonic psychometric, at-bound baseline) at a glance.
%
% Inputs:
%   tbl_trials  - long-form trial table (raw, ideally pre-demean)
%   fig_dir     - output directory
%   subj_labels - cell array of subject ID strings
%   baselines   - optional: per-subject baseline table from demeanTrialTable

    if nargin < 3 || isempty(subj_labels)
        subj_list = unique(tbl_trials.subject_id, 'stable');
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:numel(subj_list), 'UniformOutput', false);
    end
    if nargin < 4
        baselines = table();
    end
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    n_subj = numel(subj_list);

    % Layout: n_subj rows x 3 columns
    row_h_in = 1.9;
    fig_h = max(2.0 + n_subj * row_h_in, 7);
    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 14 fig_h], ...
        'PaperUnits', 'inches', 'PaperSize', [14 fig_h], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 14 fig_h]);
    tl = tiledlayout(n_subj, 3, 'Padding', 'loose', 'TileSpacing', 'compact');

    summary_rows = cell(n_subj, 1);

    blue = [0.20 0.45 0.75];
    green = [0.30 0.65 0.30];

    for s = 1:n_subj
        sid = subj_list(s);
        sm = tbl_trials.subject_id == sid;
        subj_trials = tbl_trials(sm, :);

        n_total = height(subj_trials);
        valid = subj_trials.response == 0 | subj_trials.response == 1;
        n_valid = sum(valid);
        response_rate = n_valid / n_total;
        accuracy = mean( ...
            (subj_trials.response(valid) == 1 & subj_trials.x_probe(valid) > 0) | ...
            (subj_trials.response(valid) == 0 & subj_trials.x_probe(valid) < 0));

        % --- Panel 1: trial counts per cell ---
        ax = nexttile(tl);
        cell_counts = nan(18, 1);
        cell_label_x = cell(18, 1);
        for c = 1:18
            [m, prev, curr] = conditionSubscriptsFromIndex(c);
            if m == 1
                manip_str = 'contrast';
                short_m = 'C';
            else
                manip_str = 'precision';
                short_m = 'P';
            end
            cell_counts(c) = sum(subj_trials.cond_manipulation == manip_str & ...
                                 subj_trials.cond_prev == prev & ...
                                 subj_trials.cond_curr == curr);
            cell_label_x{c} = sprintf('%s%d%d', short_m, prev, curr);
        end
        b = bar(ax, 1:18, cell_counts);
        b.FaceColor = 'flat';
        for c = 1:18
            [m, ~, ~] = conditionSubscriptsFromIndex(c);
            if m == 1
                b.CData(c, :) = blue;
            else
                b.CData(c, :) = green;
            end
        end
        xline(ax, 9.5, 'k-', 'HandleVisibility', 'off');
        set(ax, 'TickDir', 'out', 'XTick', 1:18, 'XTickLabel', cell_label_x, ...
                'XTickLabelRotation', 90, 'FontSize', 7);
        ylabel(ax, sprintf('%s\nn trials', subj_labels{s}), 'FontSize', 9, 'Interpreter', 'none');
        if s == 1
            title(ax, 'Trials per condition cell', 'FontSize', 9);
        end
        box(ax, 'off');

        % --- Panel 2: empirical psychometric ---
        ax = nexttile(tl);
        x_unique = unique(round(subj_trials.x_probe(valid)));
        p_cw_emp = nan(numel(x_unique), 1);
        n_at_x = zeros(numel(x_unique), 1);
        for xi = 1:numel(x_unique)
            mask_x = abs(subj_trials.x_probe - x_unique(xi)) < 0.5 & valid;
            n_at_x(xi) = sum(mask_x);
            if n_at_x(xi) >= 5
                p_cw_emp(xi) = mean(subj_trials.response(mask_x), 'omitnan');
            end
        end
        sz = 80 * (n_at_x / max(max(n_at_x), 1));
        sz(sz < 10) = 10;
        scatter(ax, x_unique, p_cw_emp, sz, 'filled', 'MarkerFaceColor', [0.4 0.4 0.4], ...
                'MarkerEdgeColor', 'k');
        hold(ax, 'on');
        % overlay fitted psychometric if baseline available
        if ~isempty(baselines) && ismember('subject_id', baselines.Properties.VariableNames)
            row = baselines(baselines.subject_id == sid, :);
            if ~isempty(row) && isfinite(row.mu_baseline) && isfinite(row.sigma_baseline)
                xfit = linspace(min(x_unique) - 1, max(x_unique) + 1, 200);
                p_psy = (1 - 0.25) * normcdf(xfit, row.mu_baseline, row.sigma_baseline) + 0.5 * 0.25;
                plot(ax, xfit, p_psy, 'r-', 'LineWidth', 1.2);
                xline(ax, row.mu_baseline, 'r--', 'HandleVisibility', 'off');
            end
        end
        ylim(ax, [0, 1]);
        yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');
        xline(ax, 0, 'k:', 'HandleVisibility', 'off');
        xlabel(ax, 'x_probe (deg)', 'FontSize', 8, 'Interpreter', 'none');
        ylabel(ax, 'P(CW)', 'FontSize', 8);
        if s == 1
            title(ax, 'Empirical psychometric', 'FontSize', 9);
        end
        set(ax, 'TickDir', 'out', 'FontSize', 7);
        box(ax, 'off');

        % --- Panel 3: summary text ---
        ax = nexttile(tl);
        axis(ax, 'off');

        bl_mu = NaN; bl_sigma = NaN; skip_str = '';
        if ~isempty(baselines) && ismember('subject_id', baselines.Properties.VariableNames)
            row = baselines(baselines.subject_id == sid, :);
            if ~isempty(row)
                bl_mu = row.mu_baseline;
                bl_sigma = row.sigma_baseline;
                if ismember('skipped', baselines.Properties.VariableNames) && row.skipped
                    skip_str = '  [DEMEAN SKIPPED]';
                end
            end
        end

        txt = {
            sprintf('Subject %s', subj_labels{s});
            '';
            sprintf('N trials:        %d', n_total);
            sprintf('Response rate:   %.1f%%', 100 * response_rate);
            sprintf('Overall accuracy: %.1f%%', 100 * accuracy);
            sprintf('Median cell N:   %d', median(cell_counts));
            sprintf('Min cell N:      %d', min(cell_counts));
            '';
            sprintf('Baseline mu:     %s deg', formatNumOrNaN(bl_mu, '%+.2f'));
            sprintf('Baseline sigma:  %s deg%s', formatNumOrNaN(bl_sigma, '%.2f'), skip_str)
            };
        text(ax, 0.05, 0.95, strjoin(txt, newline), 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
             'FontSize', 7.5, 'FontName', 'Courier');
        if s == 1
            title(ax, 'Summary stats', 'FontSize', 9);
        end

        summary_rows{s} = {sid, string(subj_labels{s}), n_total, n_valid, ...
                           accuracy, response_rate, median(cell_counts), min(cell_counts), ...
                           bl_mu, bl_sigma};
    end

    title(tl, 'Per-subject raw-data quality diagnostic', 'FontSize', 11);
    out_pdf = fullfile(fig_dir, 'subject_data_quality.pdf');
    try
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    catch exportErr
        warning('plotSubjectDataQuality:vectorExportFailed', ...
            'Vector export failed (%s). Falling back to raster PDF.', exportErr.message);
        exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 200);
    end
    close(fig);

    summary_table = cell2table(vertcat(summary_rows{:}), ...
        'VariableNames', {'subject_id', 'subject_label', 'n_trials', 'n_valid_responses', ...
                          'accuracy', 'response_rate', 'median_cell_n', 'min_cell_n', ...
                          'baseline_mu', 'baseline_sigma'});
    writetable(summary_table, fullfile(fig_dir, 'subject_data_quality_summary.csv'));
end

function s = formatNumOrNaN(x, fmt)
    if isnan(x); s = 'NA'; else; s = sprintf(fmt, x); end
end
