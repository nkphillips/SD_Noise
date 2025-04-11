%%% simulate_subject
close all; clear all; clc

commandwindow;

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

save_data = 1;

%% Set directories

script_dir = pwd;
functions_dir = '../functions'; addpath(functions_dir);
data_dir = '../data'; addpath(data_dir);

%% Set subjects

num_subjs = 10;

subj_IDs = cell(num_subjs,1);

for subj = 1:num_subjs
    subj_IDs{subj} = num2str(900 + subj - 1);
end

%% Simulate subjects

for subj = 1:num_subjs

    %% Set subject ID and folder

    p.subj_ID = subj_IDs{subj};

    if ~exist([data_dir '/' p.subj_ID],'dir')
        mkdir([data_dir '/' p.subj_ID])
    end

    %% Define stimuli parameters

    cond_names = {'contrast', 'filter'};
    p.num_conds = numel(cond_names);
    % Define contrasts
    p.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

    % Define orientation bandpass filter widths
    p.bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °

    if length(p.bp_filter_width) == length(p.contrast)
        p.num_levels = length(p.contrast);
    else
        disp('Condition levels do not match in length!');
    end

    % Define orientation range
    p.orientation_min = 0;
    p.orientation_max = 179;

    %% Define experiment

    p.num_blocks = 6;
    p.num_blocks_per_cond = p.num_blocks / p.num_conds;

    p.block_order = repmat(1:p.num_conds, 1, p.num_blocks_per_cond);
    p.block_order = Shuffle(p.block_order);

    first_block = nan(1, p.num_conds);
    for cond = 1:p.num_conds
        first_block(cond) = find(p.block_order == cond, 1, 'first');
    end

    p.num_trials_per_cond = 200;  

    % Generate level order, orientations, correct response
    for n_block = 1:p.num_blocks
    
        if p.block_order(n_block) == 1
            
            level_order = BalanceFactors(p.num_trials_per_cond, 1, 1:p.num_levels);
            
        elseif p.block_order(n_block) == 2
            
            level_order = BalanceFactors(p.num_trials_per_cond, 1, 1:p.num_levels);
            
        end

        if n_block == 1
            p.num_trials_per_block = length(level_order);
            p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks); % num_trials x [test_orientation, probe_orientation, cond_lvl] x num_blocks
            p.correct_response = nan(p.num_trials_per_block, p.num_blocks);
        end

        % Sample Test orientations
        test_orientation = sample_orientation(p.orientation_min, p.orientation_max, p.num_trials_per_block);

        % Storing trial events
        p.trial_events(:,:,n_block) = [test_orientation, nan(length(test_orientation),1), level_order];
        test_orientation_col = 1;
        probe_orientation_col = 2;
        level_order_col = 3;

    end

    p.num_trials = p.num_trials_per_block * p.num_blocks;

    %% Simulate responses

    mu_range = [-0.1 0.1];
    sigma_range = [2 4];
    guess_rate_range = [0.01 0.05];

    mu = mu_range(1) + (mu_range(2) - mu_range(1)) .* rand(1);
    sigma = sigma_range(1) + (sigma_range(2) - sigma_range(1)) .* rand(1);
    guess_rate = guess_rate_range(1) + (guess_rate_range(2) - guess_rate_range(1)) .* rand(1);

    p.psychometric_params = [mu, sigma, guess_rate];

    [behav_data.responses, behav_data.correct] = simulate_responses(p);

    %% Save the data

    save_filename = ['SD_Noise_S' p.subj_ID '_simulated.mat'];

    run_info.behav_data = behav_data;
    run_info.p = p;

    cd(data_dir)

    if ~exist(p.subj_ID,'dir')
        mkdir(p.subj_ID)
        disp(['/' p.subj_ID ' created.'])
    end

    if save_data
        save([p.subj_ID '/' save_filename],'run_info','-mat','-v7.3');
    end

    cd(script_dir)

end