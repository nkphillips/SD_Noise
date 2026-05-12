function out = plotCurrentOrientationErrorSummary(tbl_trials, fig_dir, subj_labels, ps)
% plotCurrentOrientationErrorSummary  Error rate as a function of current orientation.
%
% Error is defined from the 2AFC probe decision: a CW response is correct when
% x_probe > 0 and a CCW response is correct when x_probe < 0. Zero-offset trials
% are excluded because correctness is undefined for the sign rule.

    out = struct('pdf', '', 'csv', '');
    if nargin < 3 || isempty(subj_labels)
        subj_list = unique(tbl_trials.subject_id, 'stable');
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:numel(subj_list), 'UniformOutput', false);
    end
    if nargin < 4 || isempty(ps)
        ps = plotSettings();
    end
    if ~ismember('current_orientation', tbl_trials.Properties.VariableNames)
        warning('plotCurrentOrientationErrorSummary:missingCurrentOrientation', ...
            'tbl_trials lacks current_orientation; rerun trial-table construction before plotting error by current orientation.');
        return
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    edges = 0:15:180;
    centers = edges(1:end-1) + diff(edges) ./ 2;
    n_bins = numel(centers);
    n_subj = numel(subj_list);
    n_manip = 2;
    manip_names = {'contrast', 'precision'};
    manip_labels = {'Contrast', 'Precision'};
    manip_colors = [ps.colors.blue; ps.colors.green];

    err_rate = nan(n_subj, n_bins, n_manip);
    n_trials = zeros(n_subj, n_bins, n_manip);

    for s = 1:n_subj
        for m = 1:n_manip
            mask = tbl_trials.subject_id == subj_list(s) & ...
                tbl_trials.cond_manipulation == manip_names{m} & ...
                isfinite(tbl_trials.current_orientation) & ...
                isfinite(tbl_trials.x_probe) & tbl_trials.x_probe ~= 0 & ...
                (tbl_trials.response == 0 | tbl_trials.response == 1);

            ori = mod(tbl_trials.current_orientation(mask), 180);
            resp = tbl_trials.response(mask);
            xp = tbl_trials.x_probe(mask);
            is_error = (resp == 1 & xp < 0) | (resp == 0 & xp > 0);

            for b = 1:n_bins
                if b == n_bins
                    in_bin = ori >= edges(b) & ori <= edges(b + 1);
                else
                    in_bin = ori >= edges(b) & ori < edges(b + 1);
                end
                n_trials(s, b, m) = sum(in_bin);
                if n_trials(s, b, m) > 0
                    err_rate(s, b, m) = mean(is_error(in_bin), 'omitnan');
                end
            end
        end
    end

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 7], ...
        'PaperUnits', 'inches', 'PaperSize', [11 7], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 7]);
    tl = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile(tl);
    hold(ax, 'on');
    for m = 1:n_manip
        y = squeeze(err_rate(:, :, m));
        y_mean = mean(y, 1, 'omitnan');
        y_sem = local_semOmitNan(y, 1);
        errorbar(ax, centers, y_mean, y_sem, '-o', ...
            'Color', manip_colors(m, :) * 0.85, ...
            'MarkerFaceColor', manip_colors(m, :), ...
            'MarkerEdgeColor', [1 1 1], ...
            'MarkerSize', 5, ...
            'LineWidth', ps.line_width, ...
            'CapSize', 0, ...
            'DisplayName', manip_labels{m});
    end
    yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');
    xlim(ax, [0 180]);
    ylim(ax, [0 1]);
    xlabel(ax, 'Current orientation (deg)', 'FontSize', ps.axes_label_font_size);
    ylabel(ax, 'Error rate', 'FontSize', ps.axes_label_font_size);
    title(ax, 'Group mean error by current orientation', 'FontSize', ps.axes_label_font_size);
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
    local_styleAxis(ax, ps);

    ax = nexttile(tl);
    hold(ax, 'on');
    for m = 1:n_manip
        n_mean = mean(squeeze(n_trials(:, :, m)), 1, 'omitnan');
        plot(ax, centers, n_mean, '-o', ...
            'Color', manip_colors(m, :) * 0.85, ...
            'MarkerFaceColor', manip_colors(m, :), ...
            'MarkerEdgeColor', [1 1 1], ...
            'MarkerSize', 5, ...
            'LineWidth', ps.line_width, ...
            'DisplayName', manip_labels{m});
    end
    xlim(ax, [0 180]);
    xlabel(ax, 'Current orientation (deg)', 'FontSize', ps.axes_label_font_size);
    ylabel(ax, 'Mean trials per subject', 'FontSize', ps.axes_label_font_size);
    title(ax, 'Orientation-bin coverage', 'FontSize', ps.axes_label_font_size);
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
    local_styleAxis(ax, ps);

    title(tl, 'Subject data quality: error by current orientation', ...
        'FontSize', ps.axes_label_font_size + 1);
    out.pdf = fullfile(fig_dir, 'current_orientation_error_summary.pdf');
    exportgraphics(fig, out.pdf, 'ContentType', 'vector');
    close(fig);

    out.csv = fullfile(fig_dir, 'current_orientation_error_summary.csv');
    local_writeSummaryCsv(out.csv, subj_list, subj_labels, centers, err_rate, n_trials, manip_names);
end

function sem = local_semOmitNan(x, dim)
    n = sum(isfinite(x), dim);
    sem = std(x, 0, dim, 'omitnan') ./ sqrt(n);
end

function local_writeSummaryCsv(csv_path, subj_list, subj_labels, centers, err_rate, n_trials, manip_names)
    subject_id = [];
    subject_label = strings(0, 1);
    manipulation = strings(0, 1);
    current_orientation_bin_center = [];
    error_rate = [];
    n = [];

    for s = 1:numel(subj_list)
        for m = 1:numel(manip_names)
            for b = 1:numel(centers)
                subject_id(end+1, 1) = subj_list(s); %#ok<AGROW>
                subject_label(end+1, 1) = string(subj_labels{s}); %#ok<AGROW>
                manipulation(end+1, 1) = string(manip_names{m}); %#ok<AGROW>
                current_orientation_bin_center(end+1, 1) = centers(b); %#ok<AGROW>
                error_rate(end+1, 1) = err_rate(s, b, m); %#ok<AGROW>
                n(end+1, 1) = n_trials(s, b, m); %#ok<AGROW>
            end
        end
    end

    T = table(subject_id, subject_label, manipulation, current_orientation_bin_center, ...
        error_rate, n);
    writetable(T, csv_path);
end

function local_styleAxis(ax, ps)
    axis(ax, 'square');
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
        'LineWidth', ps.line_width, 'FontSize', ps.axes_tick_font_size);
end
