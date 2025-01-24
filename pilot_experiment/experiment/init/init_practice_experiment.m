%%% Initialize experiment

%% Loading screen

loading_text = 'Loading practice ...'; % standby text for loading experiment and stimuli
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Prep stimuli for drawing 

init_stimuli

init_textures

make_textures

make_patches

%% Define number and sequence of events

% Giving conditions a str and an index that can be leveraged
p.hemifields = {'L' 'R'};
p.adapt_cond_names = {'No Adapt'};
p.num_conds = numel(p.adapt_cond_names);

% The minimum number of blocks is determined by the number of adaptors and hemifields
p.min_num_blocks = numel(p.adapt_cond_names) * numel(p.hemifields);
p.num_blocks = p.min_num_blocks;

%% Condition order

% Pre-define adaptor and hemifield condition order
p.adapt_cond_order = [1 1; repmat(1:numel(p.hemifields),1 , p.num_conds)]; % [adaptation condition; hemifield]

%% Define order of test sfs

p.num_reps_per_test_sf_per_block = 2;
p.num_trials_per_block = p.num_reps_per_test_sf_per_block * p.num_test_sfs;

p.delta_sf_limit = 1.5; % octaves
p.test_sf_order = nan(p.num_blocks, p.num_trials_per_block);
p.reps_per_test_sf = nan(p.num_blocks, p.num_test_sfs);
p.sf_ratio_order = nan(p.num_blocks, p.num_trials_per_block);
p.delta_sf = nan(p.num_blocks, p.num_trials_per_block-1);

for n_block = 1:p.num_blocks

    % test_sf_order = Shuffle(repmat(p.test_sfs,1,p.num_reps_per_test_sf_per_block));
    % p.test_sf_order(n_block,:) = test_sf_order;

    [test_sf_order, reps_per_test_sf] = gen_sf_seq(p.test_sfs, p.num_trials_per_block, p.delta_sf_limit, p.num_reps_per_test_sf_per_block);
    
    p.test_sf_order(n_block,:) = test_sf_order;
    p.reps_per_test_sf(n_block,:) = reps_per_test_sf;
    
    p.sf_ratio_order(n_block,:) = test_sf_order/p.reference_sf;
    p.delta_sf(n_block,:) = log2(test_sf_order(1:end-1)./test_sf_order(2:end));

end

% Check sequences
%{
figure('Color',[1 1 1])
for n_block = 1:p.num_blocks
    subplot(1,p.num_blocks,n_block), plot(p.delta_sf(n_block,:)); 
    axis square; box off; set(gca, 'TickDir','out');
    if n_block
       xlabel('Trial');ylabel('\Delta SF');
    end
end
%}

%% Initialize behavioral data struct

init_behavioral

%% Define timing and generate frame presentation

init_practice_timing

frames = init_frames(t,p);

