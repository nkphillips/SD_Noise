function plotSubjectBaselineBias(baselines, fig_dir, subj_labels)
% plotSubjectBaselineBias  Bar plots of the per-subject baseline (mu_i, sigma_i)
% used by demeanTrialTable to demean x_probe.
%
% Saves subject_baseline_bias.pdf and subject_baseline_bias.csv into fig_dir.

    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    n_subj = height(baselines);
    if nargin < 3 || isempty(subj_labels)
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:n_subj, 'UniformOutput', false);
    end

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 4.5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 4.5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 4.5]);
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Detect skipped column (added in newer demeanTrialTable)
    has_skip = ismember('skipped', baselines.Properties.VariableNames);
    if has_skip
        skipped = baselines.skipped;
    else
        skipped = false(n_subj, 1);
    end

    % Mu (baseline bias)
    ax = nexttile(tl);
    b = bar(ax, 1:n_subj, baselines.mu_baseline);
    b.FaceColor = 'flat';
    for i = 1:n_subj
        if skipped(i)
            b.CData(i, :) = [0.85 0.30 0.30];   % red for skipped
        else
            b.CData(i, :) = [0.40 0.60 0.85];   % blue for applied
        end
    end
    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    set(ax, 'XTick', 1:n_subj, 'XTickLabel', subj_labels, 'XTickLabelRotation', 45, ...
            'TickDir', 'out');
    ylabel(ax, 'baseline mu_i (deg)');
    title(ax, 'Per-subject baseline response bias mu_i (red = skipped, mu_i = 0 applied)');
    box(ax, 'off');

    % Sigma
    ax = nexttile(tl);
    b = bar(ax, 1:n_subj, baselines.sigma_baseline);
    b.FaceColor = 'flat';
    for i = 1:n_subj
        if skipped(i)
            b.CData(i, :) = [0.85 0.30 0.30];   % red for skipped (same color as mu panel)
        else
            b.CData(i, :) = [0.85 0.55 0.40];   % orange for applied
        end
    end
    set(ax, 'XTick', 1:n_subj, 'XTickLabel', subj_labels, 'XTickLabelRotation', 45, ...
            'TickDir', 'out');
    ylabel(ax, 'baseline sigma_i (deg)');
    title(ax, 'Per-subject baseline psychometric noise sigma_i (saturation at bound triggers skip)');
    box(ax, 'off');

    title(tl, 'Per-subject baseline psychometric (no DoG, ignoring Delta-theta)');

    exportgraphics(fig, fullfile(fig_dir, 'subject_baseline_bias.pdf'), 'ContentType', 'vector');
    close(fig);

    writetable(baselines, fullfile(fig_dir, 'subject_baseline_bias.csv'));
end
