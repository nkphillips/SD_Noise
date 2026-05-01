function plotBootstrapCIDiagnostics(sd, sd_ci_cluster, sd_boot_cluster, p, plt_opts, out_dir)
% plotBootstrapCIDiagnostics
%
% Diagnostic plots for subject-cluster bootstrap DoG confidence intervals.
% These figures are intended for auditing CI shape and stability, not as
% primary manuscript figures.

if nargin < 6 || isempty(out_dir)
    out_dir = pwd;
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if isempty(sd_boot_cluster) || ~isfield(sd_boot_cluster, 'params') || isempty(sd_boot_cluster.params)
    warning('plotBootstrapCIDiagnostics:noBootstrap', 'No bootstrap parameters found; skipping CI diagnostics.');
    return
end

B = size(sd_boot_cluster.params, 1);
num_levels = size(sd.all.params_est, 1);
num_conds = size(sd.all.params_est, 3);
param_info = getParamInfo(sd, sd_boot_cluster, sd_ci_cluster);

plotParameterHistograms(sd, sd_boot_cluster, sd_ci_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds);
plotCurveSpaghetti(sd, sd_boot_cluster, sd_ci_cluster, p, plt_opts, out_dir, num_levels, num_conds);
plotCIConvergence(sd, sd_boot_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds, B);
plotFitValidityAndWidth(sd_boot_cluster, sd_ci_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds);
plotBoundaryHits(sd_boot_cluster, p, plt_opts, out_dir, num_levels, num_conds);

end

function param_info = getParamInfo(sd, sd_boot_cluster, sd_ci_cluster)
num_params = min(3, size(sd_boot_cluster.params, 5));
param_info = struct('name', {}, 'file_tag', {}, 'boot', {}, 'point', {}, ...
    'ci_lo', {}, 'ci_hi', {}, ...
    'ci_pct_lo', {}, 'ci_pct_hi', {}, 'ci_bca_lo', {}, 'ci_bca_hi', {});

% Resolve percentile + BCa structures (back-compat with old saved files)
has_pct = isfield(sd_ci_cluster, 'percentile') && isstruct(sd_ci_cluster.percentile);
has_bca = isfield(sd_ci_cluster, 'bca') && isstruct(sd_ci_cluster.bca);

for param = 1:num_params
    if param == 1
        param_info(param).name = 'Amplitude';
        param_info(param).file_tag = 'amplitude';
        param_info(param).boot = ensureBootShape(sd_boot_cluster.params(:,:,:,:,1), sd);
        param_info(param).point = squeeze(sd.all.params_est(:,:,:,1));
        param_info(param).ci_lo = squeeze(sd_ci_cluster.lo(:,:,:,1));
        param_info(param).ci_hi = squeeze(sd_ci_cluster.hi(:,:,:,1));
        if has_pct
            param_info(param).ci_pct_lo = squeeze(sd_ci_cluster.percentile.lo(:,:,:,1));
            param_info(param).ci_pct_hi = squeeze(sd_ci_cluster.percentile.hi(:,:,:,1));
        end
        if has_bca
            param_info(param).ci_bca_lo = squeeze(sd_ci_cluster.bca.lo(:,:,:,1));
            param_info(param).ci_bca_hi = squeeze(sd_ci_cluster.bca.hi(:,:,:,1));
        end
    elseif param == 2
        param_info(param).name = 'FWHM';
        param_info(param).file_tag = 'fwhm';
        param_info(param).boot = ensureBootShape(sd_boot_cluster.fwhm, sd);
        param_info(param).point = (2 * sqrt(log(2))) ./ max(squeeze(sd.all.params_est(:,:,:,2)), eps);
        param_info(param).ci_lo = sd_ci_cluster.fwhm_lo;
        param_info(param).ci_hi = sd_ci_cluster.fwhm_hi;
        if has_pct
            param_info(param).ci_pct_lo = sd_ci_cluster.percentile.fwhm_lo;
            param_info(param).ci_pct_hi = sd_ci_cluster.percentile.fwhm_hi;
        end
        if has_bca
            param_info(param).ci_bca_lo = sd_ci_cluster.bca.fwhm_lo;
            param_info(param).ci_bca_hi = sd_ci_cluster.bca.fwhm_hi;
        end
    elseif param == 3
        param_info(param).name = 'Baseline';
        param_info(param).file_tag = 'baseline';
        param_info(param).boot = ensureBootShape(sd_boot_cluster.params(:,:,:,:,3), sd);
        param_info(param).point = squeeze(sd.all.params_est(:,:,:,3));
        param_info(param).ci_lo = squeeze(sd_ci_cluster.lo(:,:,:,3));
        param_info(param).ci_hi = squeeze(sd_ci_cluster.hi(:,:,:,3));
        if has_pct
            param_info(param).ci_pct_lo = squeeze(sd_ci_cluster.percentile.lo(:,:,:,3));
            param_info(param).ci_pct_hi = squeeze(sd_ci_cluster.percentile.hi(:,:,:,3));
        end
        if has_bca
            param_info(param).ci_bca_lo = squeeze(sd_ci_cluster.bca.lo(:,:,:,3));
            param_info(param).ci_bca_hi = squeeze(sd_ci_cluster.bca.hi(:,:,:,3));
        end
    end
end
end

function boot = ensureBootShape(boot, sd)
% Preserve [B, prev, curr, cond] shape even when MATLAB loads/saves a
% singleton bootstrap dimension in a squeezed form.
target_tail = size(sd.all.params_est(:,:,:,1)); % [prev, curr, cond]
if ndims(boot) == 4 && size(boot, 2) == target_tail(1) && ...
        size(boot, 3) == target_tail(2) && size(boot, 4) == target_tail(3)
    return
end
if ndims(boot) == 3 && all(size(boot) == target_tail)
    boot = reshape(boot, [1, target_tail]);
end
end

function plotParameterHistograms(sd, sd_boot_cluster, sd_ci_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds) %#ok<INUSD>
admitted_grid = []; % only for the parametric arrays (Amplitude, Baseline) sourced from sd_boot_cluster.params
if isfield(sd_boot_cluster, 'admitted')
    admitted_grid = sd_boot_cluster.admitted;
end

for param = 1:numel(param_info)
    for cond = 1:num_conds
        fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1150]);
        tiledlayout(num_levels, num_levels, 'TileSpacing', 'compact', 'Padding', 'loose');
        axs = gobjects(num_levels, num_levels);
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                nexttile;
                axs(prev_lvl, curr_lvl) = gca;
                vals_full = squeeze(param_info(param).boot(:, prev_lvl, curr_lvl, cond));
                if ~isempty(admitted_grid)
                    mask = squeeze(admitted_grid(:, prev_lvl, curr_lvl, cond));
                    vals = vals_full(mask);
                else
                    vals = vals_full;
                end
                vals = vals(isfinite(vals));
                if ~isempty(vals)
                    histogram(vals, 'FaceColor', conditionColor(cond, plt_opts), ...
                        'FaceAlpha', 0.45, 'EdgeColor', 'none');
                    hold on;
                    drawXLine(param_info(param).point(prev_lvl, curr_lvl, cond), '-', plt_opts.colors.black, 1.2);
                    drawXLine(median(vals, 'omitnan'), '-', plt_opts.colors.red, 1);

                    % Percentile CI (black dashed)
                    if ~isempty(param_info(param).ci_pct_lo)
                        drawXLine(param_info(param).ci_pct_lo(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.black, 1);
                        drawXLine(param_info(param).ci_pct_hi(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.black, 1);
                    else
                        % Back-compat: only active CI present
                        drawXLine(param_info(param).ci_lo(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.black, 1);
                        drawXLine(param_info(param).ci_hi(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.black, 1);
                    end

                    % BCa CI (red dashed)
                    if ~isempty(param_info(param).ci_bca_lo)
                        drawXLine(param_info(param).ci_bca_lo(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.red, 1);
                        drawXLine(param_info(param).ci_bca_hi(prev_lvl, curr_lvl, cond), '--', plt_opts.colors.red, 1);
                    end
                end
                title_str = levelPairTitle(cond, prev_lvl, curr_lvl, p);
                if ~isempty(admitted_grid)
                    title_str = sprintf('%s (n=%d)', title_str, numel(vals));
                end
                title(title_str, 'FontSize', 9);
                set(gca, 'TickDir', 'out', 'Box', 'off');
            end
        end
        applySharedAxisLimits(axs, true, true);
        sgtitle(sprintf('%s Bootstrap Distributions - %s (black: percentile, red: BCa)', ...
            param_info(param).name, p.cond_names{cond}));
        saveDiagnosticFigure(fg, out_dir, sprintf('CI Histogram %s %s', param_info(param).name, p.cond_names{cond}), plt_opts);
        close(fg);
    end
end
end

function drawXLine(value, style, color, width)
% xline errors on NaN/Inf; guard so missing CIs simply skip drawing
if ~isfinite(value)
    return
end
xline(value, style, 'Color', color, 'LineWidth', width);
end

function plotCurveSpaghetti(sd, sd_boot_cluster, sd_ci_cluster, p, plt_opts, out_dir, num_levels, num_conds)
if ~isfield(sd_boot_cluster, 'curves') || ~isfield(sd_boot_cluster, 'curve_x')
    return
end

B = size(sd_boot_cluster.curves, 1);
max_curves = min(100, B);
curve_sample = unique(round(linspace(1, B, max_curves)));
curve_x = sd_boot_cluster.curve_x(:)';

for cond = 1:num_conds
    fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1150]);
    tiledlayout(num_levels, num_levels, 'TileSpacing', 'compact', 'Padding', 'loose');
    axs = gobjects(num_levels, num_levels);
    for prev_lvl = 1:num_levels
        for curr_lvl = 1:num_levels
            nexttile;
            axs(prev_lvl, curr_lvl) = gca;
            hold on;
            curves = reshape(sd_boot_cluster.curves(:, prev_lvl, curr_lvl, cond, :), B, []);
            ci_lo = squeeze(sd_ci_cluster.curve_lo(prev_lvl, curr_lvl, cond, :))';
            ci_hi = squeeze(sd_ci_cluster.curve_hi(prev_lvl, curr_lvl, cond, :))';
            curve_idx = isfinite(curve_x) & isfinite(ci_lo) & isfinite(ci_hi);
            if any(curve_idx)
                fill([curve_x(curve_idx), fliplr(curve_x(curve_idx))], ...
                    [ci_lo(curve_idx), fliplr(ci_hi(curve_idx))], ...
                    conditionColor(cond, plt_opts), 'FaceAlpha', 0.20, ...
                    'EdgeColor', 'none');
            end
            for i_curve = 1:numel(curve_sample)
                curr_curve = curves(curve_sample(i_curve), :);
                if any(isfinite(curr_curve))
                    plot(curve_x, curr_curve, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
                end
            end
            dog_params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
            if all(isfinite(dog_params))
                plot(curve_x, calcDoG(curve_x, dog_params), 'Color', plt_opts.colors.black, 'LineWidth', 1.2);
            end
            line([-90 90], [0 0], 'Color', plt_opts.colors.black, 'LineWidth', 0.5);
            title(levelPairTitle(cond, prev_lvl, curr_lvl, p), 'FontSize', 9);
            xlim([-90 90]);
            set(gca, 'TickDir', 'out', 'Box', 'off');
        end
    end
    applySharedAxisLimits(axs, true, true);
    sgtitle(sprintf('Bootstrap Curve Spaghetti - %s', p.cond_names{cond}));
    saveDiagnosticFigure(fg, out_dir, sprintf('CI Curve Spaghetti %s', p.cond_names{cond}), plt_opts);
    close(fg);
end
end

function plotCIConvergence(sd, sd_boot_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds, B) %#ok<INUSD>
if B < 2
    return
end

if B <= 25
    checkpoints = 2:B;
else
    checkpoints = unique(round(linspace(max(5, ceil(B * 0.05)), B, min(20, B))));
end

for param = 1:numel(param_info)
    for cond = 1:num_conds
        fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1150]);
        tiledlayout(num_levels, num_levels, 'TileSpacing', 'compact', 'Padding', 'loose');
        axs = gobjects(num_levels, num_levels);
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                nexttile;
                axs(prev_lvl, curr_lvl) = gca;
                vals = squeeze(param_info(param).boot(:, prev_lvl, curr_lvl, cond));
                lo = nan(size(checkpoints));
                hi = nan(size(checkpoints));
                for i_chk = 1:numel(checkpoints)
                    curr_vals = vals(1:checkpoints(i_chk));
                    curr_vals = curr_vals(isfinite(curr_vals));
                    if ~isempty(curr_vals)
                        q = prctile(curr_vals, [2.5 97.5]);
                        lo(i_chk) = q(1);
                        hi(i_chk) = q(2);
                    end
                end
                plot(checkpoints, lo, '-', 'Color', conditionColor(cond, plt_opts), 'LineWidth', 1);
                hold on;
                plot(checkpoints, hi, '-', 'Color', conditionColor(cond, plt_opts), 'LineWidth', 1);
                yline(param_info(param).point(prev_lvl, curr_lvl, cond), '--', ...
                    'Color', plt_opts.colors.black, 'LineWidth', 0.75);
                title(levelPairTitle(cond, prev_lvl, curr_lvl, p), 'FontSize', 9);
                xlabel('B');
                ylabel(param_info(param).name);
                set(gca, 'TickDir', 'out', 'Box', 'off');
            end
        end
        applySharedAxisLimits(axs, true, true);
        sgtitle(sprintf('%s CI Convergence - %s', param_info(param).name, p.cond_names{cond}));
        saveDiagnosticFigure(fg, out_dir, sprintf('CI Convergence %s %s', param_info(param).name, p.cond_names{cond}), plt_opts);
        close(fg);
    end
end
end

function plotFitValidityAndWidth(sd_boot_cluster, sd_ci_cluster, param_info, p, plt_opts, out_dir, num_levels, num_conds) %#ok<INUSD>
admitted_grid = [];
if isfield(sd_boot_cluster, 'admitted')
    admitted_grid = sd_boot_cluster.admitted;
end

for param = 1:numel(param_info)
    fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 950]);
    tiledlayout(num_conds, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
    axs = gobjects(num_conds, 3);
    for cond = 1:num_conds
        vals = param_info(param).boot(:,:,:,cond);
        if ~isempty(admitted_grid)
            admit_cond = admitted_grid(:,:,:,cond); % [B, prev, curr]
            admit_prop = squeeze(mean(admit_cond, 1));
            discard_prop = 1 - admit_prop;
        else
            admit_prop   = squeeze(mean(isfinite(vals), 1));
            discard_prop = 1 - admit_prop;
        end
        ci_width = squeeze(param_info(param).ci_hi(:,:,cond) - param_info(param).ci_lo(:,:,cond));

        nexttile;
        axs(cond, 1) = gca;
        imagesc(admit_prop, [0 1]);
        colorbar;
        axis square;
        title(sprintf('%s admitted fits', p.cond_names{cond}));
        formatHeatmapAxes(cond, p);

        nexttile;
        axs(cond, 2) = gca;
        imagesc(discard_prop, [0 1]);
        colorbar;
        axis square;
        title(sprintf('%s discarded fits', p.cond_names{cond}));
        formatHeatmapAxes(cond, p);

        nexttile;
        axs(cond, 3) = gca;
        imagesc(ci_width);
        colorbar;
        axis square;
        title(sprintf('%s CI width', p.cond_names{cond}));
        formatHeatmapAxes(cond, p);
    end
    applySharedAxisLimits(axs, true, true);
    sgtitle(sprintf('%s Bootstrap Admission and CI Width', param_info(param).name));
    saveDiagnosticFigure(fg, out_dir, sprintf('CI Validity Width %s', param_info(param).name), plt_opts);
    close(fg);
end

% Per-cell discard rate by reason (one figure across all reasons)
plotDiscardByReason(sd_boot_cluster, p, plt_opts, out_dir, num_levels, num_conds); %#ok<NASGU>
end

function plotDiscardByReason(sd_boot_cluster, p, plt_opts, out_dir, num_levels, num_conds) %#ok<INUSL>
if ~isfield(sd_boot_cluster, 'discard_reason')
    return
end
reasons = sd_boot_cluster.discard_reason; % [B, prev, curr, cond] string array
if isempty(reasons)
    return
end
B = size(reasons, 1);
if B == 0
    return
end

reason_labels = {'min_windows', 'missing_lobe', 'failed_fit', 'at_bound'};
n_reasons = numel(reason_labels);

for cond = 1:num_conds
    fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 950]);
    tiledlayout(1, n_reasons, 'TileSpacing', 'compact', 'Padding', 'loose');
    axs = gobjects(1, n_reasons);
    for r = 1:n_reasons
        rate = squeeze(sum(reasons(:, :, :, cond) == reason_labels{r}, 1)) / B;
        nexttile;
        axs(1, r) = gca;
        imagesc(rate, [0 1]);
        colorbar;
        axis square;
        title(reason_labels{r}, 'Interpreter', 'none');
        formatHeatmapAxes(cond, p);
    end
    sgtitle(sprintf('Discard rate by reason - %s', p.cond_names{cond}));
    saveDiagnosticFigure(fg, out_dir, sprintf('CI Discard By Reason %s', p.cond_names{cond}), plt_opts);
    close(fg);
end
end

function plotBoundaryHits(sd_boot_cluster, p, plt_opts, out_dir, num_levels, num_conds) %#ok<INUSL>
% plotBoundaryHits
%
% Per-condition heatmap of the at-bound rate per parameter. Highlights cells
% where fmincon converged to a parameter bound, which can spike the bootstrap
% percentile and motivate filtering or BCa.
if ~isfield(sd_boot_cluster, 'at_bound')
    return
end
at_bound = sd_boot_cluster.at_bound; % [B, prev, curr, cond, num_params]
if isempty(at_bound)
    return
end
B = size(at_bound, 1);
if B == 0
    return
end
n_params = size(at_bound, 5);

param_labels = {'Amplitude', 'Width', 'Baseline', 'Sigma'};
n_params = min(n_params, numel(param_labels));

for cond = 1:num_conds
    fg = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 950]);
    tiledlayout(1, n_params, 'TileSpacing', 'compact', 'Padding', 'loose');
    axs = gobjects(1, n_params);
    for ip = 1:n_params
        rate = squeeze(mean(at_bound(:, :, :, cond, ip), 1));
        nexttile;
        axs(1, ip) = gca;
        imagesc(rate, [0 1]);
        colorbar;
        axis square;
        title(sprintf('%s', param_labels{ip}));
        formatHeatmapAxes(cond, p);
    end
    sgtitle(sprintf('At-bound rate per parameter - %s', p.cond_names{cond}));
    saveDiagnosticFigure(fg, out_dir, sprintf('CI Boundary Hits %s', p.cond_names{cond}), plt_opts);
    close(fg);
end
end

function formatHeatmapAxes(cond, p)
set(gca, 'TickDir', 'out', 'Box', 'off');
if cond == 1
    set(gca, 'XTick', 1:3, 'XTickLabel', p.contrast, 'YTick', 1:3, 'YTickLabel', p.contrast);
else
    set(gca, 'XTick', 1:3, 'XTickLabel', p.precision, 'YTick', 1:3, 'YTickLabel', p.precision);
end
xlabel('Current level');
ylabel('Previous level');
end

function title_str = levelPairTitle(cond, prev_lvl, curr_lvl, p)
if cond == 1
    title_str = sprintf('%s -> %s', p.contrast{prev_lvl}, p.contrast{curr_lvl});
else
    title_str = sprintf('%s -> %s', p.precision{prev_lvl}, p.precision{curr_lvl});
end
end

function c = conditionColor(cond, plt_opts)
if cond == 1
    c = plt_opts.colors.blue;
else
    c = plt_opts.colors.green;
end
end

function applySharedAxisLimits(axs, share_x, share_y)
axs = axs(isgraphics(axs, 'axes'));
if isempty(axs)
    return
end

has_data = arrayfun(@(ax) ~isempty(ax.Children), axs);
source_axes = axs(has_data);
if isempty(source_axes)
    source_axes = axs;
end

if share_x
    x_lims = vertcat(source_axes.XLim);
    x_lim = [min(x_lims(:,1)), max(x_lims(:,2))];
    if all(isfinite(x_lim)) && x_lim(1) < x_lim(2)
        set(axs, 'XLim', x_lim);
    end
end

if share_y
    y_lims = vertcat(source_axes.YLim);
    y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];
    if all(isfinite(y_lim)) && y_lim(1) < y_lim(2)
        set(axs, 'YLim', y_lim);
    end
end
end

function saveDiagnosticFigure(fg, out_dir, fig_name, plt_opts)
safe_name = regexprep(fig_name, '[^\w\s-]', '');
safe_name = strrep(strtrim(safe_name), ' ', '_');
out_file = fullfile(out_dir, [safe_name '.' plt_opts.fg_type]);
set(fg, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
try
    exportgraphics(fg, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    saveas(fg, out_file);
end
end
