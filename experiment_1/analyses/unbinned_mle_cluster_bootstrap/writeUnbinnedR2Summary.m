function r2_summary = writeUnbinnedR2Summary(fig_dir, summary_table)
% writeUnbinnedR2Summary  Aggregate per-cell DoG R2 metrics for reporting.
%
% summary_table already stores one R2 per previous x current x manipulation
% cell. This helper reports the mean, SD, and SEM across those cells.

    r2_summary = table();

    if nargin < 2 || isempty(summary_table) || height(summary_table) == 0
        return
    end

    rows = {};
    if ismember('r2_delta_bins', summary_table.Properties.VariableNames)
        rows = [rows; local_metricRows(summary_table, 'r2_delta_bins', ...
            'Delta-bin Efron R2')]; %#ok<AGROW>
    end
    if ismember('r2_mcfadden', summary_table.Properties.VariableNames)
        rows = [rows; local_metricRows(summary_table, 'r2_mcfadden', ...
            'Trial-level McFadden R2')]; %#ok<AGROW>
    end

    if isempty(rows)
        return
    end

    r2_summary = cell2table(rows, 'VariableNames', ...
        {'metric', 'metric_description', 'cell_subset', 'n_cells', ...
        'mean_r2', 'sd_r2', 'sem_r2'});

    if nargin >= 1 && ~isempty(fig_dir)
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir);
        end
        writetable(r2_summary, fullfile(fig_dir, 'r2_summary.csv'));
    end
end

function rows = local_metricRows(summary_table, metric_name, metric_description)
    rows = {};
    rows = [rows; local_oneRow(summary_table.(metric_name), metric_name, ...
        metric_description, 'all_cells')]; %#ok<AGROW>

    manip_vals = cellstr(string(summary_table.cond_manipulation));
    subsets = {'contrast', 'precision'};
    for i_subset = 1:numel(subsets)
        keep = strcmpi(manip_vals, subsets{i_subset});
        rows = [rows; local_oneRow(summary_table.(metric_name)(keep), metric_name, ...
            metric_description, subsets{i_subset})]; %#ok<AGROW>
    end
end

function row = local_oneRow(vals, metric_name, metric_description, cell_subset)
    vals = vals(:);
    vals = vals(isfinite(vals));
    n_cells = numel(vals);

    if n_cells == 0
        mean_r2 = NaN;
        sd_r2 = NaN;
        sem_r2 = NaN;
    else
        mean_r2 = mean(vals);
        sd_r2 = std(vals, 0);
        sem_r2 = sd_r2 ./ sqrt(n_cells);
    end

    row = {metric_name, metric_description, cell_subset, n_cells, ...
        mean_r2, sd_r2, sem_r2};
end
