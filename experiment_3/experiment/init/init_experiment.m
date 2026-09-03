%%% Initialize experiment

%% Loading screen

if toggles.training
    loading_text = ['Loading training run ' num2str(p.run_num) ' ...'];
elseif toggles.calibration
    loading_text = ['Loading calibration run ' num2str(p.run_num) ' ...'];
else
    loading_text = ['Loading experiment run ' num2str(p.run_num) ' ...'];
end

if toggles.simulation_mode == 0
    loading_text_boundary = Screen('TextBounds', w.window, loading_text);
    loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

    Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
    Screen('Flip', w.window);
end

%% Prep stimuli for drawing

init_fixation

init_stimuli_params

if toggles.simulation_mode == 0
    stimuli = init_textures(p, dirs, w, toggles);

    make_textures

    make_patches
end

%% Define number of events
% figure out how many blocks and num trials per level per block after
% changing num levels to 2

p.feature_name = {'contrast', 'filter'};
p.num_features = numel(p.feature_name);


%% FUNCTION

init_timing
[best_N, best_blocks] = recommend_experiment_design(t, p);

%%


if toggles.training
    p.num_blocks = p.num_features; % Default = num_features (2)
elseif toggles.calibration
    p.num_blocks = 8; % Default = 8
else
    p.num_blocks = best_blocks;
    while mod(p.num_blocks, p.num_features) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_features) '. Please enter a multiple of 2: ']); end
end

p.num_blocks_per_feature = p.num_blocks / p.num_features;

p.feature_order = repmat(1:p.num_features, 1, p.num_blocks_per_feature);
p.feature_order = Shuffle(p.feature_order);

first_block = nan(1, p.num_features);
for feature = 1:p.num_features
    first_block(feature) = find(p.feature_order == feature, 1, 'first');
end




if toggles.training
    % Note that the number of levels for each condition in training is 1
    p.num_trials_per_lvl_per_block = 20; % Default = 20
elseif toggles.calibration
    p.num_trials_per_lvl_per_block = 10; % Default = 10
else %change this variable
    p.num_trials_per_lvl_per_block = best_N;
end

%% Generate level order, orientations, correct response

for n_block = 1:p.num_blocks

    level_order = BalanceFactors(p.num_trials_per_lvl_per_block, 1, 1:p.num_levels);

    if n_block == 1
        p.num_trials_per_block = length(level_order);
        p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks); % num_trials x [test_orientations, probe_orientations, feature_lvl] x num_blocks
        p.correct_response = nan(p.num_trials_per_block, p.num_blocks);
    end

    % Sample orientations
    test_orientations = sample_orientation(p.orientation_min, p.orientation_max, p.num_trials_per_block);

    probe_offsets = datasample(p.probe_offsets, length(test_orientations));
    probe_orientations = calc_probe_orientation(test_orientations, probe_offsets');

    % Storing trial events
    p.trial_events(:,:,n_block) = [test_orientations, probe_orientations, level_order];
    test_orientation_col = 1;
    probe_orientation_col = 2;
    level_order_col = 3;

end

p.num_trials = p.num_trials_per_block * p.num_blocks;

%%

t.iti_dur = t.iti_min + (t.iti_max - t.iti_min) .* rand(p.num_trials-1, 1);
t.block_dur_est = round(((t.trial_dur_est + mean(t.iti_dur(1:p.num_trials_per_block)))  * p.num_trials_per_block) / 60, 2);
t.exp_dur_est = round((p.num_trials * t.trial_dur_est + sum(t.iti_dur) + t.rest_dur * (p.num_blocks-1)) / 60, 2);

if toggles.show_debug_output
    disp(['Estimated experiment duration: ' num2str(t.exp_dur_est) ' minutes (' num2str(t.block_dur_est) ' minutes / block)']);
end

%% Check condition distribution

% get_trial_distribution(p)

%% Define timing and generate frame presentation

frames = initFrames(t,p);

%% Initialize behavioral data struct

init_behavioral
