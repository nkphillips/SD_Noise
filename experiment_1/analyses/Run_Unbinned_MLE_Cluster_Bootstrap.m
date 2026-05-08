%%% Run_Unbinned_MLE_Cluster_Bootstrap
%
% Sole entry point for the trial-level (unbinned) DoG + psychometric MLE with
% subject-cluster bootstrap. Mirrors SD_Noise_Analyses_And_Figures subject list,
% setup tag, and load_matfiles; does not run the rest of the legacy pipeline.
%
% Usage: cd to experiment_1/analyses, then run this script by name.
%
% Outputs: PDFs under unbinned_mle_cluster_bootstrap/figures/; workspace variable
% `results` (or struct `results` with fields n1, n2, ... when n_back_list has
% multiple entries).

%% Prepare workspace

close all;
clear all;
clc;

%% --- User configuration (edit here) ---

which_setup = '3329C_ASUS';

analysis_date = datestr(now, 'mm.dd.yyyy');
analysis_datetime = datestr(now, 'yyyymmdd_HHMMSS');   % time-stamped subdir for figures and CSVs

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
%          Gaussian-envelope FWHM bounds (in deg) for human readability and converted.
%   sigma: psychometric noise (deg). Lower bound 1 deg avoids lapse-sigma confound at the corner.
%   beta : global response-bias offset (deg). +-10 deg accommodates observers
%          with strong CW/CCW response biases.
p.fwhm_min_deg = 8;     % DoG FWHM lower bound (deg) -- maps to w upper bound
p.fwhm_max_deg = 120;   % DoG FWHM upper bound (deg) -- maps to w lower bound
sd_w_lb = (2 * sqrt(log(2))) / p.fwhm_max_deg;
sd_w_ub = (2 * sqrt(log(2))) / p.fwhm_min_deg;

amp_ub = 30;
amp_lb = -amp_ub;

beta_ub = 10;
beta_lb = -beta_ub;

%   columns: [A,  w, sigma, beta]; rows: [upper; lower]
p.sd_bounds = [ amp_ub, sd_w_ub,  sigma_ub,  beta_ub;
amp_lb, sd_w_lb,   sigma_lb, beta_lb];

num.subjs = length(p.subj_IDs);
num.conds = length(p.cond_names);
num.levels = 3;

% Serial separation between previous and current trial (same convention as main analysis).
n_back_list = [1 2 3];

bootstrap.B = 10000;
bootstrap.seed = 1;
bootstrap.B_grid = 100;
bootstrap.use_parallel = true;
bootstrap.ci_prctile = [2.5, 97.5];
bootstrap.ci_method = 'bca';            % 'bca' | 'percentile'
bootstrap.compute_jackknife = true;     % required for BCa
bootstrap.discard_at_bound = true;
bootstrap.discard_failed_fits = true;
bootstrap.bound_tol = 1e-4;
bootstrap.guess_rate = 0.25;            % S&S 2022 fixed lapse rate
bootstrap.compute_contrasts = true;     % auto-compute within-manipulation BCa contrasts
bootstrap.contrast_params = {'A','w','sigma','beta'};      % {'A','w','sigma','beta'}

%% --- Diagnostic toggles (extra figures and CSVs alongside the main outputs) ---

toggles.subject_influence = true;            % per-subject leverage plot + heatmap (uses jackknife)
toggles.amplitude_sigma_correlation = true;  % per-subject per-manipulation DoG fits, A vs sigma scatter (S&S 2022 Fig 1H analog)
toggles.subject_data_quality = true;         % per-subject trial counts, empirical psychometric, summary stats (always uses pre-demean trials)
toggles.demean_x_probe = true;              % subtract per-subject baseline mu_i from x_probe before SD fitting (saves diagnostic CSV+plot)
toggles.skip_at_bound_baseline = true;       % when demean is on: skip demean for subjects whose baseline sigma_i saturates at a bound

% Build fit_opts struct from p.sd_bounds; passed through to fitConditionMLE
fit_opts = struct( ...
    'lb', p.sd_bounds(2, :)', ...   % 4 x 1 column vector: lower bounds [A; w; sigma; beta]
    'ub', p.sd_bounds(1, :)' ...    % 4 x 1 column vector: upper bounds
    );

% Optional: limit parallel workers ([] uses parpool default when use_parallel is true).
p.num_workers = [];

%% Paths & dependencies

addpath('functions');
addpath('unbinned_mle_cluster_bootstrap');
addpath('loading');

ps = plotSettings();

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
        'demean_x_probe', toggles.demean_x_probe, ...
        'skip_at_bound_baseline', toggles.skip_at_bound_baseline, ...
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

    if multi_nback
        results.(sprintf('n%d', n_back)) = res_i;
    else
        results = res_i;
    end
end

elapsed_s = toc(analysis_tic);
disp(' ');
fprintf('Total analysis duration: %.1f s (%.2f min)\n', elapsed_s, elapsed_s / 60);
disp('Finished. Vector figures: unbinned_mle_cluster_bootstrap/figures/ (grids + super-SD scatter PDFs).');
