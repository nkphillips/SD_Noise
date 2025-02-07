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

p.num_blocks = 6;
while mod(p.num_blocks, p.num_conds) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_conds) '. Please enter a multiple of 2: ']); end

p.num_blocks_per_cond = p.num_blocks / p.num_conds;

p.block_order = repmat(1:p.num_conds, 1, p.num_blocks_per_cond);
p.block_order = Shuffle(p.block_order);

if p.training 
    p.num_trials_per_unique_cond = p.num_levels;  
else
    p.num_trials_per_unique_cond = 1005; % default = 105
end

p.num_trials_per_block = p.num_levels * p.num_trials_per_unique_cond/p.num_blocks_per_cond;
p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks); % num_trials x [test_orientation, probe_orientation, cond_lvl] x num_blocks 
p.correct_response = nan(p.num_trials_per_block, p.num_blocks);

for n_block = 1:p.num_blocks
   
    if p.block_order(n_block) == 1
        
        level_order = BalanceFactors(p.num_trials_per_block/p.num_levels, 1, 1:length(stimuli.contrast));
        
    elseif p.block_order(n_block) == 2
        
        level_order = BalanceFactors(p.num_trials_per_block/p.num_levels, 1, 1:length(stimuli.bp_filter_width));
        
    end
    
    test_orientation = round(stimuli.orientation_min + (stimuli.orientation_max - stimuli.orientation_min) .* rand(p.num_trials_per_block, 1));
    probe_orientation = round(stimuli.orientation_min + (stimuli.orientation_max - stimuli.orientation_min) .* rand(p.num_trials_per_block, 1));
    
    
    % Storing trial events
    p.trial_events(:,:,n_block) = [test_orientation, probe_orientation, level_order];
   
    
    % Storing correct response
    correct_response = double(test_orientation < probe_orientation);
    correct_response(correct_response == 0) = 2;
    p.correct_response(:,n_block) = correct_response;

end

p.num_trials = size(p.trial_events,1) * size(p.trial_events,3);
p.num_trials = p.num_trials_per_block * p.num_blocks;

%% Check condition distribution

% "scoreboard" for the # of level pairs
pair_count = zeros(p.num_levels, p.num_levels, p.num_conds);

% store the indices of the current trial if the current and previous trial
% match the current level pair
% note: you must use a cell array
pair_indices = cell(p.num_levels, p.num_levels, p.num_blocks_per_cond, p.num_conds);

for cond = 1:p.num_conds
    
    curr_blocks = find(p.block_order == cond);
    
    curr_trial_events = squeeze(p.trial_events(:,end,curr_blocks));
    
    % Count # of unique trial pairs (ignore 1st trial)
    for level_a = 1:p.num_levels % prev
        for level_b = 1:p.num_levels % curr
            
            % use level_a and level_b to identify the current trial pair (trial_n and trial_{n-1})
            % count the number of those pairs across each block
            for n_block = 1:length(curr_blocks)
                for n_trial = 2:size(curr_trial_events,1)
                    
                    prev_trial_lvl = curr_trial_events(n_trial-1, n_block);
                    curr_trial_lvl = curr_trial_events(n_trial, n_block);
                    
                    if prev_trial_lvl == level_a  && curr_trial_lvl == level_b
                        
                        pair_count(level_a, level_b, cond) = pair_count(level_a, level_b) + 1;
                        
                        % store the indices of the current trial if the current and previous trial
                        pair_indices{level_a, level_b, n_block, cond} = [pair_indices{level_a, level_b, n_block, cond}, n_trial];
                   
                    end
                end
            end            
        end
    end
    
end

% test for diff in number of pairs between condition

% x = pair_count(:,:,1);
% y = pair_count(:,:,2);
% 
% [h, p_value] = ttest(x(:), y(:));
% 
% figure, histogram(x(:)), hold on;  histogram(y(:))

%% Define timing and generate frame presentation

init_timing

frames = init_frames(t,p);

%% Initialize behavioral data struct

init_behavioral


