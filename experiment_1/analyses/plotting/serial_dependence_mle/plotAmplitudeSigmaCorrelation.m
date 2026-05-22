function r_summary = plotAmplitudeSigmaCorrelation(T, fig_dir)
% plotAmplitudeSigmaCorrelation  S&S 2022 Fig 1H analog:
% scatter of per-subject DoG amplitude A vs psychometric noise sigma, one point per
% (subject, manipulation) pair. At-bound fits are plotted with open markers and
% excluded from the headline correlation (with the all-points correlation also reported).
% Returns a struct with the various Pearson r and p values.

    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    blue  = [0.20 0.45 0.75];
    green = [0.30 0.65 0.30];
    colors = [blue; green];

    valid = isfinite(T.A) & isfinite(T.sigma);
    Tv = T(valid, :);

    if ismember('at_bound_any', T.Properties.VariableNames)
        good = ~Tv.at_bound_any;
    else
        good = true(height(Tv), 1);
    end

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 5], ...
        'PaperUnits', 'inches', 'PaperSize', [11 5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 5]);
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    manip_names = {'contrast', 'precision'};
    r_per_manip       = nan(1, 2);  p_per_manip       = nan(1, 2);  n_per_manip       = zeros(1, 2);
    r_per_manip_good  = nan(1, 2);  p_per_manip_good  = nan(1, 2);  n_per_manip_good  = zeros(1, 2);

    % --- Left panel: A vs sigma ---
    ax = nexttile(tl);
    hold(ax, 'on');
    h = gobjects(2, 1);
    for m = 1:2
        rows = Tv.manipulation == manip_names{m};
        sg = Tv.sigma(rows);
        Av = Tv.A(rows);
        bnd = ~good(rows);

        % plot good (filled) and at-bound (open) points
        if any(~bnd)
            h(m) = scatter(ax, sg(~bnd), Av(~bnd), 80, colors(m, :), 'filled', ...
                           'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
        if any(bnd)
            scatter(ax, sg(bnd), Av(bnd), 80, colors(m, :), ...
                    'MarkerEdgeColor', colors(m, :), 'LineWidth', 1.5);
        end

        n_per_manip(m) = numel(sg);
        n_per_manip_good(m) = sum(~bnd);

        if numel(sg) >= 3
            [r_m, p_m] = corr(sg, Av, 'rows', 'complete');
            r_per_manip(m) = r_m;
            p_per_manip(m) = p_m;
        end

        sg_g = sg(~bnd); Av_g = Av(~bnd);
        if numel(sg_g) >= 3
            [r_g, p_g] = corr(sg_g, Av_g, 'rows', 'complete');
            r_per_manip_good(m) = r_g;
            p_per_manip_good(m) = p_g;
            lin = polyfit(sg_g, Av_g, 1);
            xfit = linspace(min(sg_g), max(sg_g), 50);
            plot(ax, xfit, polyval(lin, xfit), '-', 'Color', colors(m, :), 'LineWidth', 1.8);
        end
    end

    sg_pool_g = Tv.sigma(good); A_pool_g = Tv.A(good);
    if numel(sg_pool_g) >= 3
        [r_pooled_good, p_pooled_good] = corr(sg_pool_g, A_pool_g, 'rows', 'complete');
    else
        r_pooled_good = NaN; p_pooled_good = NaN;
    end
    [r_pooled_all, p_pooled_all] = corr(Tv.sigma, Tv.A, 'rows', 'complete');

    text(ax, 0.04, 0.95, sprintf('contrast (good only, n=%d): r=%.2f, p=%.3f', ...
         n_per_manip_good(1), r_per_manip_good(1), p_per_manip_good(1)), ...
         'Units', 'normalized', 'Color', blue, 'FontWeight', 'bold');
    text(ax, 0.04, 0.88, sprintf('precision (good only, n=%d): r=%.2f, p=%.3f', ...
         n_per_manip_good(2), r_per_manip_good(2), p_per_manip_good(2)), ...
         'Units', 'normalized', 'Color', green, 'FontWeight', 'bold');
    text(ax, 0.04, 0.78, sprintf('pooled (good only, n=%d): r=%.2f, p=%.3f', ...
         numel(sg_pool_g), r_pooled_good, p_pooled_good), ...
         'Units', 'normalized', 'Color', 'k', 'FontWeight', 'bold');
    text(ax, 0.04, 0.68, sprintf('all points (n=%d): r=%.2f, p=%.3f', ...
         height(Tv), r_pooled_all, p_pooled_all), ...
         'Units', 'normalized', 'Color', [0.5 0.5 0.5], 'FontSize', 9);
    if any(~good)
        text(ax, 0.04, 0.04, 'Open markers = fits at parameter bound (excluded from headline r)', ...
             'Units', 'normalized', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
    end

    legend(ax, h(isgraphics(h)), manip_names(isgraphics(h)), 'Location', 'southeast');
    xlabel(ax, '\sigma (deg)', 'Interpreter', 'tex');
    ylabel(ax, 'A (deg)');
    title(ax, 'A vs sigma per (subject, manipulation)', 'Interpreter', 'none');
    box(ax, 'off'); set(ax, 'TickDir', 'out');

    % --- Right panel: FWHM vs sigma ---
    ax = nexttile(tl);
    hold(ax, 'on');
    valid_f = isfinite(T.fwhm) & isfinite(T.sigma);
    Tvf = T(valid_f, :);
    if ismember('at_bound_any', T.Properties.VariableNames)
        good_f = ~Tvf.at_bound_any;
    else
        good_f = true(height(Tvf), 1);
    end
    r_fwhm_per_manip_good = nan(1, 2); p_fwhm_per_manip_good = nan(1, 2);
    for m = 1:2
        rows = Tvf.manipulation == manip_names{m};
        sg = Tvf.sigma(rows);
        Fv = Tvf.fwhm(rows);
        bnd = ~good_f(rows);
        if any(~bnd)
            scatter(ax, sg(~bnd), Fv(~bnd), 80, colors(m, :), 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
        if any(bnd)
            scatter(ax, sg(bnd), Fv(bnd), 80, colors(m, :), ...
                    'MarkerEdgeColor', colors(m, :), 'LineWidth', 1.5);
        end
        sg_g = sg(~bnd); Fv_g = Fv(~bnd);
        if numel(sg_g) >= 3
            [r_g, p_g] = corr(sg_g, Fv_g, 'rows', 'complete');
            r_fwhm_per_manip_good(m) = r_g;
            p_fwhm_per_manip_good(m) = p_g;
            lin = polyfit(sg_g, Fv_g, 1);
            xfit = linspace(min(sg_g), max(sg_g), 50);
            plot(ax, xfit, polyval(lin, xfit), '-', 'Color', colors(m, :), 'LineWidth', 1.8);
        end
    end
    text(ax, 0.04, 0.95, sprintf('contrast (good): r=%.2f, p=%.3f', ...
         r_fwhm_per_manip_good(1), p_fwhm_per_manip_good(1)), ...
         'Units', 'normalized', 'Color', blue, 'FontWeight', 'bold');
    text(ax, 0.04, 0.88, sprintf('precision (good): r=%.2f, p=%.3f', ...
         r_fwhm_per_manip_good(2), p_fwhm_per_manip_good(2)), ...
         'Units', 'normalized', 'Color', green, 'FontWeight', 'bold');
    xlabel(ax, '\sigma (deg)', 'Interpreter', 'tex');
    ylabel(ax, 'FWHM (deg)');
    title(ax, 'FWHM vs sigma per (subject, manipulation)', 'Interpreter', 'none');
    box(ax, 'off'); set(ax, 'TickDir', 'out');

    title(tl, 'Per-subject DoG fits across manipulations (open markers = at-bound, excluded from headline r)', ...
          'Interpreter', 'none');
    exportgraphics(fig, fullfile(fig_dir, 'amplitude_sigma_correlation.pdf'), 'ContentType', 'vector');
    close(fig);

    writetable(T, fullfile(fig_dir, 'per_subject_per_manipulation_fits.csv'));

    r_summary = struct();
    r_summary.r_pooled_all     = r_pooled_all;
    r_summary.p_pooled_all     = p_pooled_all;
    r_summary.r_pooled_good    = r_pooled_good;
    r_summary.p_pooled_good    = p_pooled_good;
    r_summary.r_per_manip      = r_per_manip;
    r_summary.p_per_manip      = p_per_manip;
    r_summary.r_per_manip_good = r_per_manip_good;
    r_summary.p_per_manip_good = p_per_manip_good;
    r_summary.n_per_manip      = n_per_manip;
    r_summary.n_per_manip_good = n_per_manip_good;
    r_summary.r_fwhm_per_manip_good = r_fwhm_per_manip_good;
    r_summary.p_fwhm_per_manip_good = p_fwhm_per_manip_good;
    r_summary.n_at_bound_total      = sum(~good);
end
