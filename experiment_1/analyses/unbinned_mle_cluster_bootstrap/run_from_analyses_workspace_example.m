% RUN_FROM_ANALYSES_WORKSPACE_EXAMPLE
%
% Example pattern (run from experiment_1/analyses with workspace populated):
%
%   addpath('functions');
%   addpath('loading');
%   addpath('unbinned_mle_cluster_bootstrap');
%
%   % After SD_Noise-style setup: analysis_date, data_dir, p.subj_IDs, num.* , load_matfiles
%   n_back = 1;
%   tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back);
%   results = runUnbinnedMLEClusterBootstrap(tbl_trials, 'B', 2000, 'seed', 1);
%
% Reduce B for a quick smoke test, e.g. 'B', 5, 'use_parallel', false.
