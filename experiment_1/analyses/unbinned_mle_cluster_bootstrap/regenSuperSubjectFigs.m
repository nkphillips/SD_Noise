function out = regenSuperSubjectFigs(sd_noise, opts)
% regenSuperSubjectFigs  Re-render super-subject unbinned MLE figures from sd_noise.
%
% This function performs no model fitting/bootstrap. It only consumes saved
% results and writes figures/CSVs.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = local_defaults(sd_noise, opts);

    if ~exist(opts.output_root, 'dir')
        mkdir(opts.output_root);
    end

    out = struct();
    out.output_root = opts.output_root;
    out.n_back = opts.n_back_list;

    for i_ci = 1:numel(opts.ci_figure_methods)
        ci_method = opts.ci_figure_methods{i_ci};
        ci_subdir = local_ciSubdir(ci_method);

        for i_nb = 1:numel(opts.n_back_list)
            n_back = opts.n_back_list(i_nb);
            key = sprintf('n%d', n_back);
            if ~isfield(sd_noise.results, key) || isempty(sd_noise.results.(key))
                warning('regenSuperSubjectFigs:missingResult', 'Skipping %s (missing result).', key);
                continue
            end

            res = local_withCIMethod(sd_noise.results.(key), ci_method, n_back);
            local_renderSuperSubjectResult(sd_noise, key, res, opts, n_back, ci_subdir);

            folded_key = sprintf('n%d_folded_delta_theta', n_back);
            if opts.regenerate_folded_delta_theta && isfield(sd_noise.results, folded_key) && ...
                    ~isempty(sd_noise.results.(folded_key))
                res_folded = local_withCIMethod(sd_noise.results.(folded_key), ci_method, n_back);
                local_renderSuperSubjectResult(sd_noise, folded_key, res_folded, opts, n_back, ...
                    fullfile(ci_subdir, 'folded_delta_theta'));
            end
        end
    end

    if numel(opts.n_back_list) > 1
        for i_ci = 1:numel(opts.ci_figure_methods)
            ci_method = opts.ci_figure_methods{i_ci};
            ci_results = local_resultsWithCIMethod(sd_noise.results, opts.n_back_list, ci_method);
            cross_dir = fullfile(opts.output_root, 'across_n_back', local_ciSubdir(ci_method));
            plotUnbinnedAmplitudeByNBack(cross_dir, ci_results, opts.n_back_list, ...
                opts.contrast_labels, opts.precision_labels, opts.ps, opts.ci_prctile, ci_method);
            plotUnbinnedFwhmByNBack(cross_dir, ci_results, opts.n_back_list, ...
                opts.contrast_labels, opts.precision_labels, opts.ps, opts.ci_prctile, ci_method);
            plotUnbinnedSdScatterByNBack(cross_dir, ci_results, opts.n_back_list, ...
                opts.contrast_labels, opts.precision_labels, opts.ps, opts.ci_prctile, ci_method);
            plotTargetedDoGEndpointEffectsByNBack(cross_dir, ci_results, opts.n_back_list, ...
                opts.contrast_labels, opts.precision_labels, opts.ps, ...
                'show_middle_level', opts.show_middle_level_endpoint);
            trend_dir = fullfile(cross_dir, 'all_level_trends');
            trend_tests = local_collectAllLevelTrendTests(ci_results, opts.n_back_list);
            if ~isempty(trend_tests) && height(trend_tests) > 0
                if ~exist(trend_dir, 'dir')
                    mkdir(trend_dir);
                end
                writetable(trend_tests, fullfile(trend_dir, 'all_level_trend_tests.csv'));
                plotUnbinnedAllLevelTrendTestsByNBack(trend_dir, trend_tests, opts.ps, ci_method);
            end
            simple_dir = fullfile(cross_dir, 'simple_slope_trends');
            simple_tests = local_collectSimpleSlopeTrendTests(ci_results, opts.n_back_list);
            if ~isempty(simple_tests) && height(simple_tests) > 0
                if ~exist(simple_dir, 'dir')
                    mkdir(simple_dir);
                end
                writetable(simple_tests, fullfile(simple_dir, 'simple_slope_trend_tests.csv'));
                plotUnbinnedSimpleSlopeTrendTestsByNBack(simple_dir, simple_tests, opts.ps, ci_method, ...
                    opts.contrast_labels, opts.precision_labels);
                writeUnbinnedSimpleSlopeTrendHtml(simple_dir, simple_tests, ci_method);
            end
        end
    end

    try
        out.statistical_results = writeUnbinnedStatisticalResultsHtml(sd_noise, opts.output_root, opts);
    catch reportErr
        warning('regenSuperSubjectFigs:statisticalReportFailed', ...
            'Statistical results HTML report failed: %s', reportErr.message);
    end

end

function opts = local_defaults(sd_noise, opts)
    if ~isfield(opts, 'ps') || isempty(opts.ps)
        opts.ps = plotSettings();
    end
    if ~isfield(opts, 'n_back_list') || isempty(opts.n_back_list)
        opts.n_back_list = sd_noise.config.n_back_list;
    end
    if ~isfield(opts, 'cond_names') || isempty(opts.cond_names)
        opts.cond_names = sd_noise.config.p.cond_names;
    end
    if ~isfield(opts, 'contrast_labels') || isempty(opts.contrast_labels)
        opts.contrast_labels = sd_noise.config.p.contrast;
    end
    if ~isfield(opts, 'precision_labels') || isempty(opts.precision_labels)
        opts.precision_labels = sd_noise.config.p.precision;
    end
    if ~isfield(opts, 'mu_bounds') || isempty(opts.mu_bounds)
        opts.mu_bounds = [sd_noise.config.p.rb_bounds(2, 1), sd_noise.config.p.rb_bounds(1, 1)];
    end
    if ~isfield(opts, 'ci_prctile') || isempty(opts.ci_prctile)
        opts.ci_prctile = sd_noise.config.bootstrap.ci_prctile;
    end
    if ~isfield(opts, 'output_root') || isempty(opts.output_root)
        opts.output_root = fullfile('unbinned_mle_cluster_bootstrap', 'figures', ...
            sd_noise.meta.analysis_datetime);
    end
    if ~isfield(opts, 'regenerate_folded_delta_theta') || isempty(opts.regenerate_folded_delta_theta)
        opts.regenerate_folded_delta_theta = true;
    end
    if ~isfield(opts, 'show_middle_level_endpoint') || isempty(opts.show_middle_level_endpoint)
        opts.show_middle_level_endpoint = false;
    end
    if ~isfield(opts, 'ci_figure_methods') || isempty(opts.ci_figure_methods)
        opts.ci_figure_methods = {'percentile', 'bca'};
    elseif ischar(opts.ci_figure_methods) || isstring(opts.ci_figure_methods)
        opts.ci_figure_methods = cellstr(opts.ci_figure_methods);
    end
end

function results_out = local_resultsWithCIMethod(results_in, n_back_list, ci_method)
    results_out = results_in;
    for i_nb = 1:numel(n_back_list)
        key = sprintf('n%d', n_back_list(i_nb));
        if isfield(results_out, key) && ~isempty(results_out.(key))
            results_out.(key) = local_withCIMethod(results_out.(key), ci_method, n_back_list(i_nb));
        end
    end
end

function res = local_withCIMethod(res, ci_method, n_back)
    if nargin < 3 || isempty(n_back)
        n_back = NaN;
    end
    ci_method = lower(char(ci_method));
    ci_field = local_ciField(ci_method);
    if ~isfield(res, ci_field) || isempty(res.(ci_field))
        warning('regenSuperSubjectFigs:missingCI', ...
            'Result is missing %s; keeping existing active CI.', ci_field);
        return
    end

    res.ci_method = ci_method;
    res.ci_active = res.(ci_field);
    if isfield(res, 'summary_table') && ~isempty(res.summary_table)
        res.summary_table = local_summaryWithActiveCI(res.summary_table, res.ci_active, ci_method);
        res.r2_summary = writeUnbinnedR2Summary([], res.summary_table);
    end
    try
        res.summary_table_ci_diagnostics = buildUnbinnedCIDiagnostics(res);
    catch
        res.summary_table_ci_diagnostics = table();
    end

    if isfield(res, 'close_far_sigma') && isfield(res.close_far_sigma, ci_field) && ...
            ~isempty(res.close_far_sigma.(ci_field))
        res.close_far_sigma.ci_active = res.close_far_sigma.(ci_field);
        if isfield(res.close_far_sigma, 'summary_table') && ~isempty(res.close_far_sigma.summary_table)
            res.close_far_sigma.summary_table = local_closeFarSummaryWithActiveCI( ...
                res.close_far_sigma.summary_table, res.close_far_sigma.ci_active);
        end
    end

    res.contrast_table = table();
    try
        cspecs = buildStandardContrasts('params', {'A','w','sigma','beta'});
        res.contrast_table = computeUnbinnedContrasts(res, cspecs);
    catch
    end
    try
        res.all_level_trend_tests = computeUnbinnedAllLevelTrendTests(res, n_back);
    catch
        res.all_level_trend_tests = table();
    end
    try
        res.simple_slope_trend_tests = computeUnbinnedSimpleSlopeTrendTests(res, n_back);
    catch
        res.simple_slope_trend_tests = table();
    end
end

function trend_tests = local_collectAllLevelTrendTests(results, n_back_list)
    trend_tests = table();
    for i_nb = 1:numel(n_back_list)
        n_back = n_back_list(i_nb);
        key = sprintf('n%d', n_back);
        if ~isfield(results, key) || isempty(results.(key))
            continue
        end
        res = results.(key);
        if isfield(res, 'all_level_trend_tests') && ~isempty(res.all_level_trend_tests)
            tbl = res.all_level_trend_tests;
        else
            try
                tbl = computeUnbinnedAllLevelTrendTests(res, n_back);
            catch
                tbl = table();
            end
        end
        if ~isempty(tbl) && height(tbl) > 0
            trend_tests = [trend_tests; tbl]; %#ok<AGROW>
        end
    end
end

function simple_tests = local_collectSimpleSlopeTrendTests(results, n_back_list)
    simple_tests = table();
    for i_nb = 1:numel(n_back_list)
        n_back = n_back_list(i_nb);
        key = sprintf('n%d', n_back);
        if ~isfield(results, key) || isempty(results.(key))
            continue
        end
        res = results.(key);
        if isfield(res, 'simple_slope_trend_tests') && ~isempty(res.simple_slope_trend_tests)
            tbl = res.simple_slope_trend_tests;
        else
            try
                tbl = computeUnbinnedSimpleSlopeTrendTests(res, n_back);
            catch
                tbl = table();
            end
        end
        if ~isempty(tbl) && height(tbl) > 0
            simple_tests = [simple_tests; tbl]; %#ok<AGROW>
        end
    end
end

function tbl = local_summaryWithActiveCI(tbl, active, ci_method)
    tbl.ci_method = repmat(string(ci_method), height(tbl), 1);
    tbl.A_ci_lo = active.param_lo(:, 1);
    tbl.A_ci_hi = active.param_hi(:, 1);
    tbl.w_ci_lo = active.param_lo(:, 2);
    tbl.w_ci_hi = active.param_hi(:, 2);
    tbl.sigma_ci_lo = active.param_lo(:, 3);
    tbl.sigma_ci_hi = active.param_hi(:, 3);
    tbl.beta_ci_lo = active.param_lo(:, 4);
    tbl.beta_ci_hi = active.param_hi(:, 4);
    tbl.fwhm_ci_lo = active.fwhm_lo;
    tbl.fwhm_ci_hi = active.fwhm_hi;
end

function tbl = local_closeFarSummaryWithActiveCI(tbl, active)
    if isfield(active, 'lo') && isfield(active, 'hi') && size(active.lo, 2) >= 3
        lo_names = {'sigma_close_ci_lo', 'sigma_far_ci_lo', 'delta_sigma_ci_lo'};
        hi_names = {'sigma_close_ci_hi', 'sigma_far_ci_hi', 'delta_sigma_ci_hi'};
        for k = 1:numel(lo_names)
            if ismember(lo_names{k}, tbl.Properties.VariableNames)
                tbl.(lo_names{k}) = active.lo(:, k);
            end
            if ismember(hi_names{k}, tbl.Properties.VariableNames)
                tbl.(hi_names{k}) = active.hi(:, k);
            end
        end
    end
end

function ci_field = local_ciField(ci_method)
    switch lower(char(ci_method))
        case 'percentile'
            ci_field = 'ci_percentile';
        case 'bca'
            ci_field = 'ci_bca';
        otherwise
            error('regenSuperSubjectFigs:badCIMethod', ...
                'ci_figure_methods entries must be ''percentile'' or ''bca''.');
    end
end

function subdir = local_ciSubdir(ci_method)
    switch lower(char(ci_method))
        case 'percentile'
            subdir = 'percentile_ci';
        case 'bca'
            subdir = 'bca_ci';
        otherwise
            subdir = [char(ci_method) '_ci'];
    end
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

function sd = local_getDerivedSd(sd_noise, key, res)
    if isfield(res, 'summary_table') && ~isempty(res.summary_table)
        sd = deriveSdStructFromUnbinnedResult(res);
    elseif isfield(sd_noise, 'derived') && isfield(sd_noise.derived, key) && ...
            isfield(sd_noise.derived.(key), 'sd')
        sd = sd_noise.derived.(key).sd;
    else
        sd = deriveSdStructFromUnbinnedResult(res);
    end
end

function local_renderSuperSubjectResult(sd_noise, key, res, opts, n_back, result_subdir)
    result_dir = fullfile(opts.output_root, sprintf('%d_back', n_back));
    if ~isempty(result_subdir)
        result_dir = fullfile(result_dir, result_subdir);
    end

    super_dir = fullfile(result_dir, 'super_subject');
    if ~exist(super_dir, 'dir')
        mkdir(super_dir);
    end

    if isfield(res, 'summary_table') && ~isempty(res.summary_table)
        writetable(res.summary_table, fullfile(super_dir, 'summary_table.csv'));
        if isfield(res, 'summary_table_ci_diagnostics') && ~isempty(res.summary_table_ci_diagnostics)
            writetable(res.summary_table_ci_diagnostics, fullfile(super_dir, 'summary_table_ci_diagnostics.csv'));
        end
        writeUnbinnedR2Summary(super_dir, res.summary_table);
        plotUnbinnedR2CellScatter(super_dir, res.summary_table, opts.ps, n_back);
        subject_cell_fits = local_subjectCellFits(res);
        if ~isempty(subject_cell_fits) && height(subject_cell_fits) > 0
            writetable(subject_cell_fits, fullfile(super_dir, 'subject_cell_fits.csv'));
            plotUnbinnedPooledSubjectPointSummaries(super_dir, res.summary_table, ...
                subject_cell_fits, opts.contrast_labels, opts.precision_labels, opts.ps);
        end
    end
    if isfield(res, 'contrast_table') && ~isempty(res.contrast_table)
        contrast_ci_method = local_contrastCIMethod(res);
        writetable(res.contrast_table, fullfile(super_dir, sprintf('contrasts_%s.csv', contrast_ci_method)));
        writetable(res.contrast_table, fullfile(super_dir, 'contrasts.csv'));
        if strcmp(contrast_ci_method, 'bca')
            writetable(res.contrast_table, fullfile(super_dir, 'contrasts_bca.csv'));
        end
    end
    if isfield(res, 'all_level_trend_tests') && ~isempty(res.all_level_trend_tests)
        writetable(res.all_level_trend_tests, fullfile(super_dir, 'all_level_trend_tests.csv'));
    end
    if isfield(res, 'simple_slope_trend_tests') && ~isempty(res.simple_slope_trend_tests)
        writetable(res.simple_slope_trend_tests, fullfile(super_dir, 'simple_slope_trend_tests.csv'));
    end
    if isfield(res, 'close_far_sigma') && isfield(res.close_far_sigma, 'summary_table') && ...
            ~isempty(res.close_far_sigma.summary_table)
        writetable(res.close_far_sigma.summary_table, fullfile(super_dir, 'close_far_sigma_summary.csv'));
        plotCloseFarSigmaSummary(super_dir, res.close_far_sigma.summary_table, ...
            opts.contrast_labels, opts.precision_labels, opts.ps, res.ci_prctile, n_back, res.ci_method);
    end

    is_folded = ~isempty(result_subdir) && strcmp(result_subdir, 'folded_delta_theta');
    if isfield(res, 'folded_delta_theta') && ~isempty(res.folded_delta_theta)
        is_folded = logical(res.folded_delta_theta);
    end

    plotDoGMLEBootstrapFigures(opts.ps, super_dir, res.curve_boot, res.grid, ...
        res.params_boot, res.ci_prctile, ...
        'contrast_labels', opts.contrast_labels, ...
        'precision_labels', opts.precision_labels, ...
        'mu_bounds', opts.mu_bounds, ...
        'overlay', res.overlay, ...
        'admitted', res.admitted, ...
        'curve_lo', res.ci_active.curve_lo, ...
        'curve_hi', res.ci_active.curve_hi, ...
        'ci_method', res.ci_method, ...
        'bootstrap_unit', local_bootstrapUnit(res), ...
        'folded_delta_theta', is_folded);

    plotUnbinnedSdScatterSummaries(super_dir, res.summary_table, ...
        opts.contrast_labels, opts.precision_labels, opts.ps, res.ci_prctile, res.ci_method);

    try
        tests_i = computeTargetedDoGHypothesisTests(res, n_back);
        if ~isempty(tests_i) && height(tests_i) > 0
            writetable(tests_i, fullfile(super_dir, 'targeted_dog_hypothesis_tests.csv'));
            plotTargetedDoGEndpointEffects(super_dir, res, tests_i, ...
                opts.contrast_labels, opts.precision_labels, opts.ps, ...
                'n_back', n_back, ...
                'show_middle_level', opts.show_middle_level_endpoint);
        end
    catch endpointErr
        warning('regenSuperSubjectFigs:endpointEffectsFailed', ...
            'Targeted endpoint-effect figure failed for %s: %s', key, endpointErr.message);
    end

    if is_folded
        title_tag = sprintf('Folded delta-theta unbinned MLE (95%% %s CI)', local_ciLabel(res.ci_method));
    else
        title_tag = sprintf('Unbinned MLE (95%% %s CI)', local_ciLabel(res.ci_method));
    end

    sd = local_getDerivedSd(sd_noise, key, res);
    past_opts = struct( ...
        'n_back', n_back, ...
        'cond_names', {opts.cond_names}, ...
        'title_tag', title_tag, ...
        'save', true, ...
        'fig_dir', result_dir);
    plotDoGPastLevelGrid(sd, opts.contrast_labels, opts.precision_labels, past_opts);
end

function bootstrap_unit = local_bootstrapUnit(res)
    if isfield(res, 'bootstrap_unit') && ~isempty(res.bootstrap_unit)
        bootstrap_unit = res.bootstrap_unit;
    else
        bootstrap_unit = 'subject';
    end
end

function subject_cell_fits = local_subjectCellFits(res)
    if isfield(res, 'subject_cell_fits') && ~isempty(res.subject_cell_fits)
        subject_cell_fits = res.subject_cell_fits;
    else
        subject_cell_fits = table();
    end
end

function ci_method = local_contrastCIMethod(res)
    if isfield(res, 'contrast_table') && ~isempty(res.contrast_table) && ...
            ismember('ci_method', res.contrast_table.Properties.VariableNames)
        ci_method = char(res.contrast_table.ci_method(1));
    elseif isfield(res, 'ci_method') && ~isempty(res.ci_method)
        ci_method = char(res.ci_method);
    else
        ci_method = 'bca';
    end
end
