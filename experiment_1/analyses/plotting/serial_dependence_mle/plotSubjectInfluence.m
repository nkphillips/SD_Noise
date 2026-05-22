function plotSubjectInfluence(influence, fig_dir, subj_labels)
% plotSubjectInfluence  Two figures into fig_dir:
%   subject_influence_summary.pdf   -- 4 panels, one per parameter; bar chart of
%                                       median |influence| per subject; flagged
%                                       subjects highlighted in red; threshold line.
%   subject_influence_heatmap.pdf   -- 4 panels, one per parameter; n_subj x 18
%                                       cells heatmap of signed influence values.
% Also writes subject_influence_long.csv in fig_dir.

    if nargin < 3 || isempty(subj_labels)
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:size(influence.values, 1), ...
                               'UniformOutput', false);
    end
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    n_subj = size(influence.values, 1);
    num_conds = size(influence.values, 2);

    param_labels = {'A (deg)', 'w (1/deg)', '\sigma (deg)', '\beta (deg)'};
    param_files  = {'A', 'w', 'sigma', 'beta'};

    % --- Panel 1: bar chart of median |influence| per subject ---
    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 8.5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 8.5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 8.5]);
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    for k = 1:4
        ax = nexttile(tl);
        med_abs = influence.median_abs(:, k);
        flagged = influence.flagged(:, k);
        threshold = influence.flag_threshold(k);

        bars = bar(ax, 1:n_subj, med_abs);
        bars.FaceColor = 'flat';
        for i = 1:n_subj
            if flagged(i)
                bars.CData(i, :) = [0.85 0.30 0.30];
            else
                bars.CData(i, :) = [0.40 0.60 0.85];
            end
        end
        hold(ax, 'on');
        if isfinite(threshold)
            yline(ax, threshold, 'r--', sprintf('flag threshold = %.3f', threshold), ...
                  'LabelHorizontalAlignment', 'right', 'FontSize', 8);
        end
        set(ax, 'XTick', 1:n_subj, 'XTickLabel', subj_labels, 'XTickLabelRotation', 45, ...
                'TickDir', 'out');
        ylabel(ax, sprintf('median |theta_{full} - theta_{LOO,i}|'), 'Interpreter', 'tex');
        title(ax, param_labels{k}, 'Interpreter', 'tex');
        box(ax, 'off');
    end
    title(tl, sprintf(['Per-subject leverage on each parameter (median |influence| across %d cells, ', ...
                       'flagged in red if > median + 2 IQR)'], num_conds), ...
          'Interpreter', 'none');
    exportgraphics(fig, fullfile(fig_dir, 'subject_influence_summary.pdf'), 'ContentType', 'vector');
    close(fig);

    % --- Panel 2: heatmap of signed influence values (subjects x cells) ---
    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 12 9], ...
        'PaperUnits', 'inches', 'PaperSize', [12 9], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 12 9]);
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    for k = 1:4
        ax = nexttile(tl);
        v = squeeze(influence.values(:, :, k));
        m_abs = max(abs(v(:)), [], 'omitnan');
        if isnan(m_abs) || m_abs == 0; m_abs = 1; end
        imagesc(ax, v, [-m_abs, m_abs]);
        colormap(ax, parula);
        cb = colorbar(ax);
        cb.Label.String = 'theta_{full} - theta_{LOO,i}';
        cb.Label.Interpreter = 'tex';
        set(ax, 'YTick', 1:n_subj, 'YTickLabel', subj_labels);
        xline(ax, 9.5, 'k-', 'LineWidth', 1);   % manipulation boundary
        xlabel(ax, 'cell index (1-9 = contrast, 10-18 = precision)');
        title(ax, param_labels{k}, 'Interpreter', 'tex');
    end
    title(tl, 'Per-subject signed influence on each cell parameter', 'Interpreter', 'none');
    exportgraphics(fig, fullfile(fig_dir, 'subject_influence_heatmap.pdf'), 'ContentType', 'vector');
    close(fig);

    % --- CSV: long-form table ---
    rows = cell(n_subj * num_conds * 4, 1);
    idx = 0;
    for s = 1:n_subj
        for c = 1:num_conds
            for k = 1:4
                idx = idx + 1;
                rows{idx} = {subj_labels{s}, c, param_files{k}, ...
                             influence.values(s, c, k), abs(influence.values(s, c, k))};
            end
        end
    end
    T = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'subject_id', 'cell', 'param', 'influence', 'influence_abs'});
    writetable(T, fullfile(fig_dir, 'subject_influence_long.csv'));
end
