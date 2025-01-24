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

% Giving conditions a str and an index that can be leveraged
p.adapt_cond_names = {};
p.num_conds = numel(p.cond_names);

% The minimum number of blocks is determined by the number of adaptors and hemifields
p.min_num_blocks = numel(p.adapt_cond_names) * numel(p.hemifields);
p.num_blocks = p.min_num_blocks;


%% Define trial events


%% Initialize behavioral data struct

init_behavioral

%% Define timing and generate frame presentation

init_timing

frames = init_frames(t,p);

