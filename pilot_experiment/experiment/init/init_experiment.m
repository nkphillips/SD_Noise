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

%% Define timing and generate frame presentation

init_timing

frames = init_frames(t,p);

p.trials_per_block = round(t.block_dur/t.trial_dur_est);

%% Define number and sequence of events

p.block_cond_names = {'contrast', 'filter'};
p.num_block_conds = numel(p.block_cond_names);

p.num_blocks = 6;
while mod(p.num_blocks, p.num_block_conds) ~= 0, p.num_blocks = input(['Error! Number of blocks must be a multiple of ' num2str(p.num_block_conds) '. Please enter a multiple of 2: ']); end

p.num_blocks_per_cond = p.num_blocks/p.num_block_conds;
p.block_order = repmat(1:p.num_block_conds,1,p.num_blocks_per_cond);
p.block_order = Shuffle(p.block_order);

if p.training 
    p.num_trials_per_cond = 2;  
else
    p.num_trials_per_cond = 30;
end

for n_block = 1:p.num_blocks
   
    if p.block_order(n_block) == 1
        
        trial_order = BalanceFactors(p.num_trials_per_cond/p.num_blocks_per_cond, 0, 1:length(stimuli.contrast));
        
    elseif p.block_order(n_block) == 2
        
    end
    p.trial_events(:,n_block) = trial_order;

end

% Define block order contrast 


%% Initialize behavioral data struct

init_behavioral


