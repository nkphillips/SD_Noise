function out = regenSubjectFigs(sd_noise, opts)
% regenSubjectFigs  Re-render subject diagnostics from saved sd_noise fields.
%
% This function performs no model fitting/bootstrap. It skips diagnostics whose
% required saved inputs are absent.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = local_defaults(sd_noise, opts);

    if ~exist(opts.output_root, 'dir')
        mkdir(opts.output_root);
    end

    out = struct();
    out.output_root = opts.output_root;
    out.shared_subject_dir = fullfile(opts.output_root, 'subject');

    first_key = sprintf('n%d', opts.n_back_list(1));
    if isfield(sd_noise.results, first_key)
        first_res = sd_noise.results.(first_key);
    else
        first_res = [];
    end

    if ~isempty(first_res) && isfield(first_res, 'baseline_bias') && ~isempty(first_res.baseline_bias)
        try
            plotSubjectBaselineBias(first_res.baseline_bias, ...
                fullfile(out.shared_subject_dir, 'baseline_bias'), opts.subj_labels);
        catch plotErr
            warning('regenSubjectFigs:baselineBiasFailed', ...
                'Baseline-bias diagnostic failed: %s', plotErr.message);
        end
    end

    if isfield(sd_noise, 'trials') && isfield(sd_noise.trials, first_key) && ...
            isfield(sd_noise.trials.(first_key), 'raw') && ~isempty(sd_noise.trials.(first_key).raw)
        baselines = table();
        if ~isempty(first_res) && isfield(first_res, 'baseline_bias')
            baselines = first_res.baseline_bias;
        end
        dq_dir = fullfile(out.shared_subject_dir, 'data_quality');
        try
            plotSubjectDataQuality(sd_noise.trials.(first_key).raw, ...
                dq_dir, opts.subj_labels, baselines);
            plotCurrentOrientationErrorSummary(sd_noise.trials.(first_key).raw, ...
                dq_dir, opts.subj_labels, opts.ps);
        catch plotErr
            warning('regenSubjectFigs:dataQualityFailed', ...
                'Subject data-quality diagnostic failed: %s', plotErr.message);
        end
        try
            if ~isempty(first_res) && isfield(first_res, 'raw_performance') && ...
                    isstruct(first_res.raw_performance) && isfield(first_res.raw_performance, 'subject_overall')
                raw_performance = first_res.raw_performance;
            else
                raw_performance = computeRawPerformanceSummary(sd_noise.trials.(first_key).raw, opts.subj_labels);
            end
            plotRawPerformanceSummary(raw_performance, dq_dir, opts.ps);
        catch perfErr
            warning('regenSubjectFigs:rawPerformanceFailed', ...
                'Raw performance diagnostic failed: %s', perfErr.message);
        end
    else
        warning('regenSubjectFigs:missingRawTrials', ...
            'Skipping shared raw-performance figures because sd_noise.trials.%s.raw is missing.', first_key);
    end

    for i_nb = 1:numel(opts.n_back_list)
        n_back = opts.n_back_list(i_nb);
        key = sprintf('n%d', n_back);
        if ~isfield(sd_noise.results, key) || isempty(sd_noise.results.(key))
            warning('regenSubjectFigs:missingResult', 'Skipping %s (missing result).', key);
            continue
        end

        res = sd_noise.results.(key);
        subj_dir = fullfile(opts.output_root, sprintf('%d_back', n_back), 'subject');
        local_renderSubjectResult(res, subj_dir, opts.subj_labels, key, false);

        folded_key = sprintf('n%d_folded_delta_theta', n_back);
        if opts.regenerate_folded_delta_theta && isfield(sd_noise.results, folded_key) && ...
                ~isempty(sd_noise.results.(folded_key))
            res_folded = sd_noise.results.(folded_key);
            folded_subj_dir = fullfile(opts.output_root, sprintf('%d_back', n_back), ...
                'folded_delta_theta', 'subject');
            local_renderSubjectResult(res_folded, folded_subj_dir, opts.subj_labels, folded_key, true);
        end
    end

end

function opts = local_defaults(sd_noise, opts)
    if ~isfield(opts, 'n_back_list') || isempty(opts.n_back_list)
        opts.n_back_list = sd_noise.config.n_back_list;
    end
    if ~isfield(opts, 'subj_labels') || isempty(opts.subj_labels)
        opts.subj_labels = sd_noise.config.p.subj_IDs;
    end
    if ~isfield(opts, 'output_root') || isempty(opts.output_root)
        opts.output_root = fullfile('unbinned_mle_cluster_bootstrap', 'figures', ...
            sd_noise.meta.analysis_datetime);
    end
    if ~isfield(opts, 'regenerate_folded_delta_theta') || isempty(opts.regenerate_folded_delta_theta)
        opts.regenerate_folded_delta_theta = true;
    end
    if ~isfield(opts, 'ps') || isempty(opts.ps)
        opts.ps = plotSettings();
    end
end

function local_renderSubjectResult(res, subj_dir, subj_labels, key, render_raw_diagnostics)
    if ~exist(subj_dir, 'dir')
        mkdir(subj_dir);
    end

    if render_raw_diagnostics && isfield(res, 'baseline_bias') && ~isempty(res.baseline_bias)
        try
            plotSubjectBaselineBias(res.baseline_bias, ...
                fullfile(subj_dir, 'baseline_bias'), subj_labels);
        catch plotErr
            warning('regenSubjectFigs:baselineBiasFailed', ...
                'Baseline-bias diagnostic failed for %s: %s', key, plotErr.message);
        end
    end

    if render_raw_diagnostics && isfield(res, 'tbl_trials_raw') && ~isempty(res.tbl_trials_raw)
        baselines = table();
        if isfield(res, 'baseline_bias')
            baselines = res.baseline_bias;
        end
        dq_dir = fullfile(subj_dir, 'data_quality');
        try
            plotSubjectDataQuality(res.tbl_trials_raw, ...
                dq_dir, subj_labels, baselines);
            plotCurrentOrientationErrorSummary(res.tbl_trials_raw, ...
                dq_dir, subj_labels, plotSettings());
        catch plotErr
            warning('regenSubjectFigs:dataQualityFailed', ...
                'Subject data-quality diagnostic failed for %s: %s', key, plotErr.message);
        end
        try
            if isfield(res, 'raw_performance') && isstruct(res.raw_performance) && ...
                    isfield(res.raw_performance, 'subject_overall')
                raw_performance = res.raw_performance;
            else
                raw_performance = computeRawPerformanceSummary(res.tbl_trials_raw, subj_labels);
            end
            plotRawPerformanceSummary(raw_performance, dq_dir, plotSettings());
        catch perfErr
            warning('regenSubjectFigs:rawPerformanceFailed', ...
                'Raw performance diagnostic failed for %s: %s', key, perfErr.message);
        end
    end

    if isfield(res, 'subject_influence') && isstruct(res.subject_influence) && ...
            isfield(res.subject_influence, 'values') && ~isempty(res.subject_influence.values)
        try
            plotSubjectInfluence(res.subject_influence, ...
                fullfile(subj_dir, 'subject_influence'), subj_labels);
        catch plotErr
            warning('regenSubjectFigs:subjectInfluenceFailed', ...
                'Subject-influence diagnostic failed for %s: %s', key, plotErr.message);
        end
    end

    if isfield(res, 'per_subject_per_manip_fits') && ~isempty(res.per_subject_per_manip_fits)
        try
            plotAmplitudeSigmaCorrelation(res.per_subject_per_manip_fits, ...
                fullfile(subj_dir, 'amplitude_sigma_correlation'));
        catch plotErr
            warning('regenSubjectFigs:amplitudeSigmaFailed', ...
                'Amplitude-sigma diagnostic failed for %s: %s', key, plotErr.message);
        end
    end
end
