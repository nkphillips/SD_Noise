%%% init_behavioral

%% Pre-allocate fields

behav_data.response = nan(p.num_trials_per_block,p.num_blocks);
behav_data.correct = nan(p.num_trials_per_block, p.num_blocks);

behav_data.performance = nan(1, p.num_blocks);
