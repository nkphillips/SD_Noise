
%%% simulate_staircase

close all; clear all; clc

commandwindow;

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.disp_on = 1;
p.which_setup = 3; % 0 = MacBook, 1 = 3329B, 2 = 3329C_ASUS, 3 = S32D850
save_staircase_data = 1;

%% Set directories

p.subj_ID = '999';

dirs.script_dir = pwd;
dirs.functions_dir = '../functions'; addpath(dirs.functions_dir);
dirs.data_dir = '../data'; addpath(dirs.data_dir);

if ~exist([dirs.data_dir '/' p.subj_ID],'dir')
    mkdir([dirs.data_dir '/' p.subj_ID])
end

%% Define stimuli

cond_names = {'contrast', 'filter'};
p.num_conds = numel(cond_names);
% Define contrasts
stimuli.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

% Define orientation bandpass filter widths
stimuli.bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °

if length(stimuli.bp_filter_width) == length(stimuli.contrast)
    p.num_levels = length(stimuli.contrast);
else
    disp('Condition levels do not match in length!');
end

% Define orientation range
stimuli.orientation_min = 0;
stimuli.orientation_max = 179;

%% Define psychometric function parameters
% For each condition (contrast and filter width, 3 levels each)

p.mu = [0 0 0; 0 0 0]; 
p.sigma = [3 3 3; 3 3 3]; % Controls steepness of psychometric function
p.guess_rate = [0.02 0.02 0.02; 0.02 0.02 0.02]; % Small probability of guessing

%% Enter staircase

p.keypress_numbers = 1:2;

init_staircase

% Loop through conditions and levels
n_block = 0;
for curr_cond = 1:p.num_conds

    disp(['Condition ' num2str(cond_names{curr_cond})])

    for curr_lvl = 1:p.num_levels

        disp(['Level ' num2str(curr_lvl)])

        n_block = n_block + 1;

        for n_trial = 1:length(staircases.trial_order(:, curr_lvl, curr_cond))

            %% Get current staircase info

            % Get current test orientation
            curr_test_orient = staircases.test_orientation(n_trial,curr_lvl,curr_cond);           
            
            % Get current staircase
            curr_sc = staircases.trial_order(n_trial, curr_lvl, curr_cond);

            % Skip current trial if the max reversals have already been reached for the current staircase
            if staircases.num_reversals(curr_sc, curr_lvl, curr_cond) == staircases.max_reversals
                % continue
            end

            % Get the current staircase trial number
            curr_sc_trial = find(n_trial == staircases.trial_indices(curr_sc, :, curr_lvl, curr_cond));

            if p.disp_on
                disp(['Staircase ' num2str(curr_sc) ', trial ' num2str(curr_sc_trial)])
            end

            % Probe offset
            curr_probe_offset = staircases.probe_offsets(curr_sc,curr_sc_trial,curr_lvl,curr_cond);
            curr_probe_orientation = calc_probe_orientation(curr_test_orient,curr_probe_offset);
            orient_diff = curr_probe_orientation - curr_test_orient;
            is_CW = orient_diff > 0;    

            % Get current contrast and filter width
            if curr_cond == 1

                curr_contrast = curr_lvl;
                curr_filter_width = 1;

            elseif curr_cond == 2

                curr_contrast = 1;
                curr_filter_width = curr_lvl;

            end


            if p.disp_on
                disp(['Test Orientation: ' num2str(curr_test_orient) '°'])
                disp(['Probe Orientation: ' num2str(curr_probe_orientation) '°'])
                disp(['Orientation difference: ' num2str(orient_diff) '°'])
                if is_CW
                    disp('Correct response: Right')
                else
                    disp('Correct response: Left')
                end
            end

            %% Simulate response

            mu = p.mu(curr_cond, curr_lvl); 
            sigma = p.sigma(curr_cond, curr_lvl);
            guess_rate = p.guess_rate(curr_cond, curr_lvl); 
            
            p_CW = calc_pCW(orient_diff, mu, sigma, guess_rate);
            response = rand() < p_CW;
            correct_response = response == is_CW;
                        
            staircases.responses(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = correct_response;

            update_staircase

            if p.disp_on
                disp(' ')
            end

        end

    end
end

%% Calculate final probe offsets for each condition & level

staircases.final_probe_offsets = nan(p.num_levels, p.num_conds);

for cond = 1:p.num_conds
    for lvl = 1:p.num_levels

        % Pre-allocate vector to store the mean probe offset for a staircase
        reversal_probe_offsets = nan(staircases.num_staircases_per_cond, 1);

        % Get the contrast levels at all the reversals from both staircases
        for n_sc = 1:staircases.num_staircases_per_cond

            % Get the number of reversals
            curr_num_reversals = nnz(staircases.reversal_indices(n_sc,:,lvl,cond) == 1);

            % Get the contrast levels at the reversals
            if curr_num_reversals > 0
                tmp_reversal_probe_offsets = staircases.probe_offsets(n_sc, staircases.reversal_indices(n_sc,:,lvl,cond) == 1, lvl,cond);
            end

            % Ignore the first reversals with respect to num_reversals_to_consider (counting back from the last reversal)
            tmp_reversal_probe_offsets = tmp_reversal_probe_offsets(end-staircases.num_reversals_to_consider+1:end);

            % Average all the reversals across a staircase
            reversal_probe_offsets(n_sc) = mean(tmp_reversal_probe_offsets, 'omitnan');

        end

        % Average the two levels across staircases
        staircases.final_probe_offsets(lvl,cond) = mean(reversal_probe_offsets);

    end
end

%% Save staircase data

staircases.contrast = stimuli.contrast;
staircases.bp_filter_width = stimuli.bp_filter_width;

if save_staircase_data
    save_filename = ['staircase_data_S' p.subj_ID '_' t.the_date '_' t.the_time '.mat'];
    save([dirs.data_dir '/' p.subj_ID '/' save_filename], 'staircases');
    disp(['Staircase data saved to ' dirs.data_dir '/' p.subj_ID '/' save_filename]);
end