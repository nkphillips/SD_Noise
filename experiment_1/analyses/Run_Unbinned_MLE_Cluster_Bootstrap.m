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

p.subj_IDs = {'001', '002', '003', '004', '006', '007', '008', '009', '010', '011', '013', '015'};
p.cond_names = {'Contrast', 'Precision'};

% Labels & y-scale limits for figures (match SD_Noise_Analyses_And_Figures.m Super Subj Bias with DoG *).
p.contrast = {'90%', '50%', '25%'};
p.precision = {'2°', '40°', '80°'};
p.rb_bounds = [20, 90;   % upper bounds
    -20, 0.1]; % lower bounds

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

analysis_tic = tic;

if numel(n_back_list) == 1

    n_back = n_back_list(1);

    disp(' ');
    disp(['=== Unbinned MLE cluster bootstrap | n_back = ' num2str(n_back) ' ===']);

    tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
    disp(['tbl_trials height: ' num2str(height(tbl_trials))]);

    nv = { ...
        'B', bootstrap.B, ...
        'seed', bootstrap.seed, ...
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
        'contrast_labels', p.contrast, ...
        'precision_labels', p.precision, ...
        'mu_bounds', [p.rb_bounds(2, 1), p.rb_bounds(1, 1)], ...
        'fig_subdir', sprintf('%d_back', n_back) ...
        };

    if ~isempty(p.num_workers)
        nv = [nv, {'num_workers', p.num_workers}]; %#ok<AGROW>
    end

    results = runUnbinnedMLEClusterBootstrap(tbl_trials, nv{:});

else

    results = struct();
    for i_nb = 1:numel(n_back_list)

        n_back = n_back_list(i_nb);

        disp(' ');
        disp(['=== Unbinned MLE cluster bootstrap | n_back = ' num2str(n_back) ' ===']);

        tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
        disp(['tbl_trials height: ' num2str(height(tbl_trials))]);

        nv = { ...
            'B', bootstrap.B, ...
            'seed', bootstrap.seed + i_nb - 1, ...
            'B_grid', bootstrap.B_grid, ...
            'use_parallel', bootstrap.use_parallel, ...
            'ci_prctile', bootstrap.ci_prctile, ...
            'ci_method', bootstrap.ci_method, ...
            'compute_jackknife', bootstrap.compute_jackknife, ...
            'discard_at_bound', bootstrap.discard_at_bound, ...
            'discard_failed_fits', bootstrap.discard_failed_fits, ...
            'bound_tol', bootstrap.bound_tol, ...
            'guess_rate', bootstrap.guess_rate, ...
            'contrast_labels', p.contrast, ...
            'precision_labels', p.precision, ...
            'mu_bounds', [p.rb_bounds(2, 1), p.rb_bounds(1, 1)], ...
            'fig_subdir', sprintf('%d_back', n_back) ...
            };

        if ~isempty(p.num_workers)
            nv = [nv, {'num_workers', p.num_workers}]; %#ok<AGROW>
        end

        results.(sprintf('n%d', n_back)) = runUnbinnedMLEClusterBootstrap(tbl_trials, nv{:});
    end

end

elapsed_s = toc(analysis_tic);
disp(' ');
fprintf('Total analysis duration: %.1f s (%.2f min)\n', elapsed_s, elapsed_s / 60);
disp('Finished. Vector figures: unbinned_mle_cluster_bootstrap/figures/ (grids + super-SD scatter PDFs).');
