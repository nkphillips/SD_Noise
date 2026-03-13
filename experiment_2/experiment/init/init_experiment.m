%%% Initialize experiment

%% Loading screen

if p.training
    loading_text = ['Loading training run ' num2str(p.run_num) ' ...'];
elseif p.calibration
    loading_text = ['Loading calibration run ' num2str(p.run_num) ' ...'];
else
    loading_text = ['Loading experiment run ' num2str(p.run_num) ' ...'];
end

loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Prep stimuli for drawing

init_fixation

init_stimuli_params

stimuli = init_textures(p, dirs, w);

make_textures

make_patches

%% Define number and sequence of events

p.feature_name = {'contrast', 'filter'};
p.num_features = numel(p.feature_name);

if p.training
    p.num_blocks = p.num_features;
elseif p.calibration
    p.num_blocks = 8;
else
    p.num_blocks = 6;
    while mod(p.num_blocks, p.num_features) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_features) '. Please enter a multiple of 2: ']); end
end

p.num_blocks_per_feature = p.num_blocks / p.num_features;

p.feature_order = repmat(1:p.num_features, 1, p.num_blocks_per_feature);
p.feature_order = Shuffle(p.feature_order);

first_block = nan(1, p.num_features);
for feature = 1:p.num_features
    first_block(feature) = find(p.feature_order == feature, 1, 'first');
end

if p.training
    % Note that the number of levels for each condition in training is 1
    p.num_trials_per_lvl_per_block = 20; % Default = 20
elseif p.calibration
    p.num_trials_per_lvl_per_block = 1; % default = 16; Based on 90m session recommendation
else
    p.num_trials_per_lvl_per_block = 40; % Default = 40
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

%% Check condition distribution

% get_trial_distribution(p)

%% Define timing and generate frame presentation

init_timing
frames = initFrames(t,p);

%% Initialize behavioral data struct

init_behavioral


