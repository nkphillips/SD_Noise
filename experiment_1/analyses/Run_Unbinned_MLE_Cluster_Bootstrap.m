%%% Run_Unbinned_MLE_Cluster_Bootstrap
%
% Sole entry point for the trial-level (unbinned) DoG + psychometric MLE with
% subject-cluster bootstrap. Mirrors SD_Noise_Analyses_And_Figures subject list,
% setup tag, and load_matfiles; does not run the rest of the legacy pipeline.
%
% Usage: cd to experiment_1/analyses, then run this script by name.
%
% Outputs: PDFs under unbinned_mle_cluster_bootstrap/figures/; per-n-back
% estimate snapshots plus a canonical sd_noise struct under
% unbinned_mle_cluster_bootstrap/estimates/<analysis_datetime>/; workspace
% variables `results` and `sd_noise`.

%% Prepare workspace

close all;
clear all;
clc;

%% --- User configuration (edit here) ---

which_setup = '3329C_ASUS';

analysis_date = datestr(now, 'mm.dd.yyyy');
analysis_datetime = datestr(now, 'yyyymmdd_HHMMSS');   % time-stamped subdir for figures and CSVs

%% Paths & dependencies

addpath('functions');
addpath('unbinned_mle_cluster_bootstrap');
addpath('plotting');
addpath('loading');

% Figure-only regeneration from a saved sd_noise object. When true, this
% script skips data loading, MLE, bootstrap, jackknife, and estimate saving.
regenerate_figures_only = false;
regenerate_analysis_datetime = '';   % leave empty to use newest saved sd_noise file
regenerate_sd_noise_file = '';
regenerate_overwrite_original_figures = false;
regenerate_super_subject = true;
regenerate_subject = true;
regenerate_folded_delta_theta = true;
show_middle_level_endpoint = false;  % faded, unconnected L2 reference dot in endpoint-effect figures

p.subj_IDs = {'001', '002', '003', '004', '006', '007', '008', '009', '010', '011', '013', '015'};
p.cond_names = {'Contrast', 'Precision'};

% Labels & y-scale limits for figures (match SD_Noise_Analyses_And_Figures.m Super Subj Bias with DoG *).
p.contrast = {'90%', '50%', '25%'};
p.precision = {'2°', '40°', '80°'};

% Response-bias bounds (used as the y-axis range for figure ribbons; not for SD model fitting)
mu_ub = 20;
mu_lb = -mu_ub;

sigma_ub = 90;
sigma_lb = 1;

%   columns: [mu, sigma]; rows: [upper; lower]
p.rb_bounds = [mu_ub, sigma_ub;
    mu_lb,sigma_lb];

% --- DoG-model bounds for the unbinned MLE (used for fmincon in fitConditionMLE) ---
% Order: [A (deg); w (1/deg); sigma (deg); beta (deg)]
%   A    : DoG peak amplitude. Allow attractive AND repulsive serial dependence.
%   w    : DoG width (1/deg) under the S&S 2022 parameterization. Specified via
%          numerically estimated DoG-lobe FWHM bounds (in deg) for human readability.
%   sigma: psychometric noise (deg). Lower bound 1 deg avoids lapse-sigma confound at the corner.
%   beta : global response-bias offset (deg). +-10 deg accommodates observers
%          with strong CW/CCW response biases.
fit_defaults = unbinnedMLEFitDefaults();
p.fwhm_min_deg = fit_defaults.fwhm_min_deg;
p.fwhm_max_deg = fit_defaults.fwhm_max_deg;

%   columns: [A,  w, sigma, beta]; rows: [upper; lower]
p.sd_bounds = [fit_defaults.ub';
    fit_defaults.lb'];

num.subjs = length(p.subj_IDs);
num.conds = length(p.cond_names);
num.levels = 3;

% Serial separation between previous and current trial (same convention as main analysis).
n_back_list = [1 2 3];

%% bootstrap settings

bootstrap.B = 10000;
bootstrap.seed = 1;
bootstrap.B_grid = 100;
bootstrap.use_parallel = true;
bootstrap.unit = 'subject';              % 'subject' primary inference | 'trial' descriptive diagnostic
bootstrap.ci_prctile = [2.5, 97.5];
bootstrap.ci_method = 'bca';            % 'bca' | 'percentile'
bootstrap.compute_jackknife = true;     % required for BCa
bootstrap.discard_at_bound = true;
bootstrap.discard_failed_fits = true;
bootstrap.bound_tol = 1e-4;
bootstrap.guess_rate = 0.25;            % S&S 2022 fixed lapse rate
bootstrap.compute_contrasts = true;     % auto-compute within-manipulation BCa contrasts
bootstrap.contrast_params = {'A','w','sigma','beta'};      % {'A','w','sigma','beta'}

if strcmpi(bootstrap.unit, 'trial')
    bootstrap.ci_method = 'percentile';
    bootstrap.compute_jackknife = false;
end

%% --- Diagnostic toggles (extra figures and CSVs alongside the main outputs) ---

toggles.subject_influence = true;            % per-subject leverage plot + heatmap (uses jackknife)
toggles.amplitude_sigma_correlation = true;  % per-subject per-manipulation DoG fits, A vs sigma scatter (S&S 2022 Fig 1H analog)
toggles.subject_data_quality = true;         % per-subject trial counts, empirical psychometric, summary stats (always uses pre-demean trials)
toggles.demean_x_probe = true;              % subtract per-subject baseline mu_i from x_probe before SD fitting (saves diagnostic CSV+plot)
toggles.skip_at_bound_baseline = true;       % when demean is on: skip demean for subjects whose baseline sigma_i saturates at a bound
toggles.fold_delta_theta_fit = false;        % also fit a folded odd-symmetric model: -Delta trials reflected to +Delta
toggles.close_far_sigma = true;              % fixed-DoG diagnostic: sigma for |Delta-theta| < threshold vs far trials
toggles.subject_cell_points = true;          % pooled point estimates over faint per-subject per-cell fits

close_far_sigma.threshold_deg = 30;
close_far_sigma.min_trials = 10;

% Fit options passed through to fitConditionMLE.
fit_opts = fit_defaults.fit_opts;

% Optional: limit parallel workers ([] uses parpool default when use_parallel is true).
p.num_workers = [];

ps = plotSettings();

if regenerate_figures_only
    analysis_tic = tic;

    sd_noise_file = regenerate_sd_noise_file;
    if isempty(sd_noise_file)
        estimates_root = fullfile('unbinned_mle_cluster_bootstrap', 'estimates');
        if isempty(regenerate_analysis_datetime)
            files = dir(fullfile(estimates_root, '*', 'SD_Noise_Unbinned_sd_noise_*.mat'));
        else
            files = dir(fullfile(estimates_root, regenerate_analysis_datetime, ...
                'SD_Noise_Unbinned_sd_noise_*.mat'));
        end
        if isempty(files)
            error('Run_Unbinned_MLE_Cluster_Bootstrap:noSavedAnalysis', ...
                'No saved sd_noise files found under %s.', estimates_root);
        end
        [~, ord] = sort([files.datenum], 'descend');
        files = files(ord);
        sd_noise_file = fullfile(files(1).folder, files(1).name);
    end

    S = load(sd_noise_file, 'sd_noise');
    if ~isfield(S, 'sd_noise')
        error('Run_Unbinned_MLE_Cluster_Bootstrap:missingSdNoise', ...
            'File does not contain variable sd_noise: %s', sd_noise_file);
    end
    sd_noise = S.sd_noise;
    fprintf('Loaded sd_noise: %s\n', sd_noise_file);

    if regenerate_overwrite_original_figures
        output_root = sd_noise.paths.figure_root;
    else
        regen_datetime = datestr(now, 'yyyymmdd_HHMMSS');
        output_root = fullfile(sd_noise.paths.figure_root, ['regen_' regen_datetime]);
    end

    if ~exist(output_root, 'dir')
        mkdir(output_root);
    end

    regen_opts = struct( ...
        'output_root', output_root, ...
        'ps', ps, ...
        'n_back_list', sd_noise.config.n_back_list, ...
        'regenerate_folded_delta_theta', regenerate_folded_delta_theta, ...
        'show_middle_level_endpoint', show_middle_level_endpoint);

    if regenerate_super_subject
        regenSuperSubjectFigs(sd_noise, regen_opts);
    end

    if regenerate_subject
        regenSubjectFigs(sd_noise, regen_opts);
    end

    elapsed_s = toc(analysis_tic);
    fprintf('Regenerated figures under: %s\n', output_root);
    fprintf('Figure regeneration duration: %.1f s (%.2f min)\n', elapsed_s, elapsed_s / 60);
    return
end

init_paths

%% Load experiment runs

load_matfiles

%% Run unbinned pipeline per n_back
% Single source of truth: one nv block, one loop. When n_back_list is a single
% value, `results` is the bare results struct from runUnbinnedMLEClusterBootstrap.
% When n_back_list has multiple values, `results` is a struct with fields n1, n2,
% ..., each holding the corresponding per-n_back results. Adding new options
% requires editing exactly one nv block.

analysis_tic = tic;
multi_nback = numel(n_back_list) > 1;

if multi_nback
    results = struct();
end

sd_noise = struct();
sd_noise.meta = struct( ...
    'pipeline', 'unbinned_mle_cluster_bootstrap', ...
    'analysis_date', analysis_date, ...
    'analysis_datetime', analysis_datetime, ...
    'created_at', datestr(now), ...
    'which_setup', which_setup);
sd_noise.config = struct( ...
    'p', p, ...
    'num', num, ...
    'bootstrap', bootstrap, ...
    'toggles', toggles, ...
    'close_far_sigma', close_far_sigma, ...
    'fit_opts', fit_opts, ...
    'n_back_list', n_back_list, ...
    'which_setup', which_setup);
sd_noise.paths = struct( ...
    'analysis_dir', pwd, ...
    'figure_root', fullfile('unbinned_mle_cluster_bootstrap', 'figures', analysis_datetime), ...
    'estimate_root', fullfile('unbinned_mle_cluster_bootstrap', 'estimates', analysis_datetime));
sd_noise.results = struct();
sd_noise.trials = struct();
sd_noise.derived = struct();

for i_nb = 1:numel(n_back_list)

    n_back = n_back_list(i_nb);
    seed_b = bootstrap.seed + (i_nb - 1) * multi_nback;   % offset only when looping multiple n_back

    disp(' ');
    disp(['=== Unbinned MLE cluster bootstrap | n_back = ' num2str(n_back) ' ===']);

    tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
    disp(['tbl_trials height: ' num2str(height(tbl_trials))]);

    nv = { ...
        'B', bootstrap.B, ...
        'seed', seed_b, ...
        'B_grid', bootstrap.B_grid, ...
        'use_parallel', bootstrap.use_parallel, ...
        'bootstrap_unit', bootstrap.unit, ...
        'ci_prctile', bootstrap.ci_prctile, ...
        'ci_method', bootstrap.ci_method, ...
        'compute_jackknife', bootstrap.compute_jackknife, ...
        'discard_at_bound', bootstrap.discard_at_bound, ...
        'discard_failed_fits', bootstrap.discard_failed_fits, ...
        'bound_tol', bootstrap.bound_tol, ...
        'guess_rate', bootstrap.guess_rate, ...
        'compute_contrasts', bootstrap.compute_contrasts, ...
        'contrast_params', bootstrap.contrast_params, ...
        'fit_opts', fit_opts, ...
        'subject_influence', toggles.subject_influence, ...
        'amplitude_sigma_correlation', toggles.amplitude_sigma_correlation, ...
        'subject_data_quality', toggles.subject_data_quality, ...
        'make_figures', true, ...
        'write_raw_subject_diagnostics', i_nb == 1, ...
        'raw_subject_fig_subdir', analysis_datetime, ...
        'n_back', n_back, ...
        'demean_x_probe', toggles.demean_x_probe, ...
        'skip_at_bound_baseline', toggles.skip_at_bound_baseline, ...
        'close_far_sigma', toggles.close_far_sigma, ...
        'close_far_sigma_threshold_deg', close_far_sigma.threshold_deg, ...
        'close_far_sigma_min_trials', close_far_sigma.min_trials, ...
        'subject_cell_points', toggles.subject_cell_points, ...
        'subj_labels', p.subj_IDs, ...
        'contrast_labels', p.contrast, ...
        'precision_labels', p.precision, ...
        'mu_bounds', [p.rb_bounds(2, 1), p.rb_bounds(1, 1)], ...
        'fig_subdir', fullfile(analysis_datetime, sprintf('%d_back', n_back)) ...
        };

    if ~isempty(p.num_workers)
        nv = [nv, {'num_workers', p.num_workers}]; %#ok<AGROW>
    end

    res_i = runUnbinnedMLEClusterBootstrap(tbl_trials, nv{:});

    result_key = sprintf('n%d', n_back);
    est_dir = fullfile('unbinned_mle_cluster_bootstrap', 'estimates', ...
        analysis_datetime, sprintf('%d_back', n_back));
    if ~exist(est_dir, 'dir')
        mkdir(est_dir);
    end
    est_fname = sprintf('SD_Noise_Unbinned_Estimates_%s_%s.mat', ...
        analysis_datetime, bootstrap.ci_method);
    sd = struct('all', struct('params_est', packSummaryTableToSdParamsEst(res_i.summary_table)));
    summary_table = res_i.summary_table;
    close_far_sigma_result = res_i.close_far_sigma;
    meta = struct( ...
        'pipeline', 'unbinned_mle_cluster_bootstrap', ...
        'analysis_datetime', analysis_datetime, ...
        'which_setup', which_setup, ...
        'n_back', n_back, ...
        'ci_method', bootstrap.ci_method, ...
        'bootstrap_unit', res_i.bootstrap_unit, ...
        'guess_rate', bootstrap.guess_rate, ...
        'n_subj', res_i.n_subj, ...
        'B', res_i.B, ...
        'delta_theta_centers', res_i.empirical_delta_centers);
    save(fullfile(est_dir, est_fname), 'sd', 'meta', 'summary_table', 'close_far_sigma_result', '-v7.3');
    fprintf('Saved unbinned estimates: %s\n', fullfile(est_dir, est_fname));

    sd_noise.results.(result_key) = res_i;
    sd_noise.trials.(result_key) = struct( ...
        'raw', res_i.tbl_trials_raw, ...
        'demeaned', res_i.tbl_trials);
    sd_noise.derived.(result_key) = struct( ...
        'sd', sd, ...
        'summary_table', summary_table, ...
        'meta', meta);

    if multi_nback
        results.(result_key) = res_i;
    else
        results = res_i;
    end

    if toggles.fold_delta_theta_fit
        disp(' ');
        disp(['=== Folded delta-theta unbinned MLE | n_back = ' num2str(n_back) ' ===']);

        [tbl_trials_folded, fold_info] = foldTrialTableDeltaTheta(tbl_trials);
        disp(['folded tbl_trials height: ' num2str(height(tbl_trials_folded))]);
        disp(['negative-delta trials flipped: ' num2str(fold_info.n_flipped)]);

        folded_fig_subdir = fullfile(analysis_datetime, sprintf('%d_back', n_back), 'folded_delta_theta');
        nv_folded = nv;
        nv_names = nv_folded(1:2:end);
        idx_fig_name = find(strcmp(nv_names, 'fig_subdir'), 1);
        idx_fig_subdir = (idx_fig_name - 1) * 2 + 1;
        nv_folded{idx_fig_subdir + 1} = folded_fig_subdir;
        idx_raw_name = find(strcmp(nv_names, 'raw_subject_fig_subdir'), 1);
        idx_raw_subdir = (idx_raw_name - 1) * 2 + 1;
        nv_folded{idx_raw_subdir + 1} = folded_fig_subdir;
        nv_folded = [nv_folded, {'folded_delta_theta', true}]; %#ok<AGROW>

        res_folded = runUnbinnedMLEClusterBootstrap(tbl_trials_folded, nv_folded{:});
        res_folded.fold_info = fold_info;

        result_key_folded = sprintf('n%d_folded_delta_theta', n_back);
        est_dir_folded = fullfile('unbinned_mle_cluster_bootstrap', 'estimates', ...
            analysis_datetime, sprintf('%d_back', n_back), 'folded_delta_theta');
        if ~exist(est_dir_folded, 'dir')
            mkdir(est_dir_folded);
        end
        est_fname_folded = sprintf('SD_Noise_Unbinned_Estimates_%s_%s_folded_delta_theta.mat', ...
            analysis_datetime, bootstrap.ci_method);
        sd_folded = struct('all', struct('params_est', packSummaryTableToSdParamsEst(res_folded.summary_table)));
        summary_table_folded = res_folded.summary_table;
        close_far_sigma_result_folded = res_folded.close_far_sigma;
        meta_folded = meta;
        meta_folded.fold_delta_theta = true;
        meta_folded.fold_info = fold_info;
        save(fullfile(est_dir_folded, est_fname_folded), ...
            'sd_folded', 'meta_folded', 'summary_table_folded', 'close_far_sigma_result_folded', '-v7.3');
        fprintf('Saved folded unbinned estimates: %s\n', fullfile(est_dir_folded, est_fname_folded));

        sd_noise.results.(result_key_folded) = res_folded;
        sd_noise.trials.(result_key_folded) = struct( ...
            'folded', res_folded.tbl_trials, ...
            'fold_source_raw', tbl_trials);
        sd_noise.derived.(result_key_folded) = struct( ...
            'sd', sd_folded, ...
            'summary_table', summary_table_folded, ...
            'meta', meta_folded);

        if multi_nback
            results.(result_key_folded) = res_folded;
        else
            results_folded_delta_theta = res_folded; %#ok<NASGU>
        end
    end
end

sd_noise.meta.completed_at = datestr(now);

regen_opts = struct( ...
    'output_root', sd_noise.paths.figure_root, ...
    'ps', ps, ...
    'n_back_list', sd_noise.config.n_back_list, ...
    'regenerate_folded_delta_theta', toggles.fold_delta_theta_fit, ...
    'show_middle_level_endpoint', show_middle_level_endpoint, ...
    'ci_figure_methods', {{'percentile', 'bca'}});
regenSuperSubjectFigs(sd_noise, regen_opts);
regenSubjectFigs(sd_noise, regen_opts);

sd_noise_path = fullfile(sd_noise.paths.estimate_root, ...
    sprintf('SD_Noise_Unbinned_sd_noise_%s_%s.mat', analysis_datetime, bootstrap.ci_method));
sd_noise.paths.sd_noise_file = sd_noise_path;
if ~exist(sd_noise.paths.estimate_root, 'dir')
    mkdir(sd_noise.paths.estimate_root);
end
save(sd_noise_path, 'sd_noise', '-v7.3');
fprintf('Saved canonical sd_noise object: %s\n', sd_noise_path);

elapsed_s = toc(analysis_tic);
disp(' ');
fprintf('Total analysis duration: %.1f s (%.2f min)\n', elapsed_s, elapsed_s / 60);
disp('Finished. Figures under unbinned_mle_cluster_bootstrap/figures/; estimates under unbinned_mle_cluster_bootstrap/estimates/ (timestamped).');
