%%% simulate_experiment

clear all; close all; clc;

%% Set directories

script_dir = pwd;
functions_dir =  '../functions'; addpath(functions_dir);

%% Toggles

p.disp_on = 1;

%% Define parameters

num_subjects = 10:5:30;
num_trials_per_cond = 40:40:200; % (1 to 5 sessions)

num_experiments = length(num_subjects) * length(num_trials_per_cond);

% Define experiment parameters (every combination of subjects and trials) 
[subject_counts, trials_per_cond_counts] = BalanceFactors(1, 0, num_subjects, num_trials_per_cond);
experiment_parameters = [subject_counts, trials_per_cond_counts];

% Define psychometric function parameters
mu_range = [-4 4]; 
sigma_range = [3 10]; 
guess_rate_range = [0.25 0.25]; 

% Define serial dependence parameter range 
% (increasing amplitude with noise levels)
% (this needs to respect each noise type and lvl pair as well)
amplitude_range = [2 4; 4 6; 6 8];
width_range = [5 10];
noise_range = [0.5 1.5];

%% Simulate experiment

experiments = struct();

for n_exp = 1:num_experiments

    % Initialize structures for current experiment
    experiments(n_exp).p = struct();
    experiments(n_exp).subj_data = struct();

    % Pull number of subjects and trials per condition for current experiment
    curr_num_subjects = experiment_parameters(n_exp, 1);
    curr_num_trials_per_cond = experiment_parameters(n_exp, 2);

    for subj = 1:curr_num_subjects

        % Define contrasts
        experiments(n_exp).p(subj).contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

        % Define orientation bandpass filter widths
        experiments(n_exp).p(subj).bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °

        % Check that the number of levels between stimulus contrast and bp filter widths match
        if length(experiments(n_exp).p(subj).bp_filter_width) == length(experiments(n_exp).p(subj).contrast)
            experiments(n_exp).p(subj).num_levels = length(experiments(n_exp).p(subj).contrast);
        else
            disp('Condition levels do not match in length!');
        end

        % Define possible orientations
        experiments(n_exp).p(subj).orientation_min = 0;
        experiments(n_exp).p(subj).orientation_max = 179;
        experiments(n_exp).p(subj).probe_offset_range = [1 15];

        % Define conditions
        experiments(n_exp).p(subj).cond_names = {'contrast', 'filter'};
        experiments(n_exp).p(subj).num_conds = numel(experiments(n_exp).p(subj).cond_names);

        % Define blocks
        experiments(n_exp).p(subj).num_blocks = 6;
        experiments(n_exp).p(subj).num_blocks_per_cond = experiments(n_exp).p(subj).num_blocks / experiments(n_exp).p(subj).num_conds;
        experiments(n_exp).p(subj).block_order = repmat(1:experiments(n_exp).p(subj).num_conds, 1, experiments(n_exp).p(subj).num_blocks_per_cond);
        experiments(n_exp).p(subj).block_order = Shuffle(experiments(n_exp).p(subj).block_order);

        % Define trial events
        experiments(n_exp).p(subj).num_trials_per_cond = curr_num_trials_per_cond;
        for n_block = 1:experiments(n_exp).p(subj).num_blocks

            if experiments(n_exp).p(subj).block_order(n_block) == 1
        
                level_order = BalanceFactors(curr_num_trials_per_cond, 1, 1:experiments(n_exp).p(subj).num_levels);
                
            elseif experiments(n_exp).p(subj).block_order(n_block) == 2
                
                level_order = BalanceFactors(curr_num_trials_per_cond, 1, 1:experiments(n_exp).p(subj).num_levels);
                
            end
        
            if n_block == 1
                experiments(n_exp).p(subj).num_trials_per_block = length(level_order);
                experiments(n_exp).p(subj).trial_events = nan(experiments(n_exp).p(subj).num_trials_per_block, 3, experiments(n_exp).p(subj).num_blocks); % num_trials x [test_orientation, probe_orientation, cond_lvl] x num_blocks
                experiments(n_exp).p(subj).correct_response = nan(experiments(n_exp).p(subj).num_trials_per_block, experiments(n_exp).p(subj).num_blocks);
            end

            test_orientation = sample_orientation(experiments(n_exp).p(subj).orientation_min, experiments(n_exp).p(subj).orientation_max, length(level_order));
            probe_offset = experiments(n_exp).p(subj).probe_offset_range(1) + (experiments(n_exp).p(subj).probe_offset_range(2) - experiments(n_exp).p(subj).probe_offset_range(1)) .* rand(length(level_order), 1);
            probe_orientation = calc_probe_orientation(test_orientation, probe_offset);

            % Store trial events
            experiments(n_exp).p(subj).trial_events(:,:,n_block) = [test_orientation, probe_orientation, level_order];

            % Define subject characteristics
            mu = mu_range(1) + (mu_range(2) - mu_range(1)) .* rand(1);
            sigma = sigma_range(1) + (sigma_range(2) - sigma_range(1)) .* rand(1);
            guess_rate = guess_rate_range(1) + (guess_rate_range(2) - guess_rate_range(1)) .* rand(1);
            experiments(n_exp).p(subj).psychometric_params = [mu, sigma, guess_rate];
            
            for i = 1:experiments(n_exp).p(subj).num_levels
                amp = amplitude_range(i,1) + (amplitude_range(i,2) - amplitude_range(i,1)) .* rand(1);
                width = width_range(1) + (width_range(2) - width_range(1)) .* rand(1);
                noise = noise_range(1) + (noise_range(2) - noise_range(1)) .* rand(1);
                experiments(n_exp).p(subj).DoG_params(i,:) = [amp width noise];
            end
    
        end

        % Simulate responses
        [experiments(n_exp).subj_data(subj).responses, experiments(n_exp).subj_data(subj).correct] = simulate_responses( experiments(n_exp).p(subj) );

        if subj == 1 && n_exp <= length(num_trials_per_cond)
            experiments(n_exp).counts(subj) = get_trial_distribution( experiments(n_exp).p(subj) );
        end

        % Summarize responses
        for cond = 1:experiments(n_exp).p(subj).num_conds
            for lvl = 1:experiments(n_exp).p(subj).num_levels

                trial_indx = experiments(n_exp).p(subj).trial_events(:,3) == lvl;
                responses = experiments(n_exp).subj_data(subj).responses(trial_indx);
                mean_response = mean(responses);

                experiments(n_exp).subj_data(subj).mean_responses(lvl, cond) = mean_response;

            end
        end

    end

    % Summarize group-level responses

    

end

%% Visualize example group response

% n_exp = find(experiment_parameters(:,1) == 10 & experiment_parameters(:,2) == p.num_levels*35);
% curr_num_subjects = experiment_parameters(n_exp,1);
% curr_num_trials_per_cond = experiment_parameters(n_exp,2);

% for cond = 1:p.num_conds

%     figure('Color', 'w','Name', ['Example group-level response for ' p.cond_names{cond}]);


% end

