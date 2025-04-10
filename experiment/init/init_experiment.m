%%% Initialize experiment

%% Loading screen

loading_text = ['Loading experiment run ' num2str(p.run_num) ' ...']; % standby text for loading experiment and stimuli
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Prep stimuli for drawing 

init_fixation

init_stimuli_params

init_textures

make_textures

make_patches

%% Define number and sequence of events

p.cond_names = {'contrast', 'filter'};
p.num_conds = numel(p.cond_names);

if p.training
    p.num_blocks = p.num_conds;
else
    p.num_blocks = 6;
    while mod(p.num_blocks, p.num_conds) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_conds) '. Please enter a multiple of 2: ']); end
end

p.num_blocks_per_cond = p.num_blocks / p.num_conds;

p.block_order = repmat(1:p.num_conds, 1, p.num_blocks_per_cond);
p.block_order = Shuffle(p.block_order);

first_block = nan(1, p.num_conds);
for cond = 1:p.num_conds
    first_block(cond) = find(p.block_order == cond, 1, 'first');
end

if p.training 
    % Note that the number of levels for each condition in training is 1
    p.num_trials_per_cond = 40;  
else
    p.num_trials_per_cond = 10; % default = 40
end

%% Generate level order, orientations, correct response

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

%% Check condition distribution

% get_trial_distribution(p)

%% Define timing and generate frame presentation

init_timing

frames = init_frames(t,p);

%% Initialize behavioral data struct

init_behavioral


