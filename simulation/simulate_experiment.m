%%% simulate_experiment

clear all; close all; clc;

%% Set directories

script_dir = pwd;
functions_dir =  '../functions'; addpath(functions_dir);

%% Toggles

p.disp_on = 1;

%% Define parameters

% Define contrasts
stimuli.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

% Define orientation bandpass filter widths
stimuli.bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °

% Check that the number of levels between stimulus contrast and bp filter widths match
if length(stimuli.bp_filter_width) == length(stimuli.contrast)
    p.num_levels = length(stimuli.contrast);
else
    disp('Condition levels do not match in length!');
end

% Define possible orientations
stimuli.orientation_min = 0;
stimuli.orientation_max = 179;
probe_offset_range = [1 15];

% Define conditions
p.cond_names = {'contrast', 'filter'};
p.num_conds = numel(p.cond_names);

num_subjects = 10:5:50;
num_trials_per_cond = p.num_levels*35:p.num_levels:p.num_levels*100;

num_experiments = length(num_subjects) * length(num_trials_per_cond);

% Define experiment parameters (every combination of subjects and trials) 
[subject_counts, trials_per_cond_counts] = BalanceFactors(1, 0, num_subjects, num_trials_per_cond);
experiment_parameters = [subject_counts, trials_per_cond_counts];

% Define psychometric function parameters
mu_range = [-2 2]; %
sigma_range = [3 10]; % Controls steepness of psychometric function
guess_rate_range = [0.01 0.25]; % 

% Define serial dependence parameter range (increasing amplitude with noise levels)
amplitude_range = [2 4; 4 6; 6 8];
width_range = [5 10];
noise_range = [0.5 1];

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

        % Define trial events
        level_order = BalanceFactors(curr_num_trials_per_cond/p.num_levels, 1, 1:p.num_levels);
        test_orientation = round(stimuli.orientation_min + (stimuli.orientation_max - stimuli.orientation_min) .* rand(curr_num_trials_per_cond, 1));
        probe_offset = probe_offset_range(1) + (probe_offset_range(2) - probe_offset_range(1)) .* rand(curr_num_trials_per_cond, 1);
        probe_orientation = calc_probe_orientation(test_orientation, probe_offset);

        % Store trial events
        experiments(n_exp).p(subj).trial_events = [test_orientation, probe_orientation, level_order];
        
        % Define subject characteristics
        mu = mu_range(1) + (mu_range(2) - mu_range(1)) .* rand(1);
        sigma = sigma_range(1) + (sigma_range(2) - sigma_range(1)) .* rand(1);
        guess_rate = guess_rate_range(1) + (guess_rate_range(2) - guess_rate_range(1)) .* rand(1);
        experiments(n_exp).p(subj).psychometric_params = [mu, sigma, guess_rate];
        
        for i = 1:p.num_levels
            amp = amplitude_range(i,1) + (amplitude_range(i,2) - amplitude_range(i,1)) .* rand(1);
            width = width_range(1) + (width_range(2) - width_range(1)) .* rand(1);
            noise = noise_range(1) + (noise_range(2) - noise_range(1)) .* rand(1);
            experiments(n_exp).p(subj).DoG_params(i,:) = [amp width noise];
        end
        
        % Simulate responses
        [experiments(n_exp).subj_data(subj).responses, ] = simulate_responses( experiments(n_exp).p(subj) );

        % Summarize responses
        for cond = 1:p.num_conds
            for lvl = 1:p.num_levels

                % Get responses for current condition and level 
                trial_indx = experiments(n_exp).p(subj).trial_events(:,3) == lvl;
                responses = experiments(n_exp).subj_data(subj).responses(trial_indx);

                % Calculate mean response for current condition and level
                mean_response = mean(responses);

                % Store mean response for current condition and level
                experiments(n_exp).subj_data(subj).mean_responses(lvl, cond) = mean_response;

            end
        end

    end

    % Summarize group-level responses

    

end



%% Visualize example group response

n_exp = find(experiment_parameters(:,1) == 10 & experiment_parameters(:,2) == p.num_levels*35);
curr_num_subjects = experiment_parameters(n_exp,1);
curr_num_trials_per_cond = experiment_parameters(n_exp,2);

for cond = 1:p.num_conds

    figure('Color', 'w','Name', ['Example group-level response for ' p.cond_names{cond}]);


end
