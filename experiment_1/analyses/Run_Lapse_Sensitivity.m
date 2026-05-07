%%% Run_Lapse_Sensitivity
%
% Sensitivity analysis: sweeps the fixed lapse rate lambda over a grid and
% reports per-cell point estimates of (A, w, FWHM, sigma, beta) at each
% lambda. NO bootstrap is performed -- the goal is to verify that
% directional conclusions about the manipulation contrasts do not depend on
% the choice of lambda. Results land under
% unbinned_mle_cluster_bootstrap/lapse_sensitivity/figures/<n_back>_back/.

%% Prepare workspace
close all; clear all; clc;

%% --- User configuration ---

which_setup = '3329C_ASUS';
analysis_date = datestr(now, 'mm.dd.yyyy');

p.subj_IDs = {'001', '002', '003', '004', '006', '007', '008', '009', '010', '011', '013', '015'};
p.cond_names = {'Contrast', 'Precision'};
p.contrast = {'90%', '50%', '25%'};
p.precision = {'2°', '40°', '80°'};
p.rb_bounds = [20, 90; -20, 1];

num.subjs = length(p.subj_IDs);
num.conds = length(p.cond_names);
num.levels = 3;

n_back_list = [1 2 3];

lapse_grid = [0.02, 0.05, 0.10, 0.25];
use_parallel = true;
p.num_workers = [];

%% Paths & dependencies
addpath('functions');
addpath('unbinned_mle_cluster_bootstrap');
addpath('unbinned_mle_cluster_bootstrap/lapse_sensitivity');
addpath('loading');

ps = plotSettings();
init_paths

%% Load runs
load_matfiles

%% Sweep
analysis_tic = tic;

if numel(n_back_list) == 1

    n_back = n_back_list(1);
    disp(' ');
    disp(['=== Lapse sensitivity sweep | n_back = ' num2str(n_back) ' ===']);

    tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
    disp(['tbl_trials height: ' num2str(height(tbl_trials))]);

    nv = { ...
        'lambdas', lapse_grid, ...
        'use_parallel', use_parallel, ...
        'fig_subdir', sprintf('%d_back', n_back), ...
        'contrast_labels', p.contrast, ...
        'precision_labels', p.precision ...
        };
    if ~isempty(p.num_workers)
        nv = [nv, {'num_workers', p.num_workers}]; %#ok<AGROW>
    end

    results = runLapseSensitivity(tbl_trials, nv{:});

else

    results = struct();
    for i_nb = 1:numel(n_back_list)
        n_back = n_back_list(i_nb);
        disp(' ');
        disp(['=== Lapse sensitivity sweep | n_back = ' num2str(n_back) ' ===']);

        tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
        disp(['tbl_trials height: ' num2str(height(tbl_trials))]);

        nv = { ...
            'lambdas', lapse_grid, ...
            'use_parallel', use_parallel, ...
            'fig_subdir', sprintf('%d_back', n_back), ...
            'contrast_labels', p.contrast, ...
            'precision_labels', p.precision ...
            };
        if ~isempty(p.num_workers)
            nv = [nv, {'num_workers', p.num_workers}]; %#ok<AGROW>
        end

        results.(sprintf('n%d', n_back)) = runLapseSensitivity(tbl_trials, nv{:});
    end

end

elapsed_s = toc(analysis_tic);
disp(' ');
fprintf('Total sensitivity sweep duration: %.1f s (%.2f min)\n', elapsed_s, elapsed_s / 60);
disp('Finished. Vector figures: unbinned_mle_cluster_bootstrap/lapse_sensitivity/figures/.');
