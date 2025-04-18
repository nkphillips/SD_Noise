% Serial dependence + noise

%{


============================
Written by Luis D. Ramirez & Natalya Phillips
lur003@ucsd.edu
UCSD
Created 10.24.2024

%}

%%

if ~p.simulate_response, wait_for_trigger; end
t.exp_start_time = GetSecs;

%% Block loop

for n_block = 1:p.num_blocks

    if p.disp_on, disp(['Block ' num2str(n_block)]); end

    % Pull the current block condition
    curr_cond = p.block_order(n_block);

    if p.use_staircase

        % Update probe offsets if it's not the first block of the current condition
        if n_block > first_block(curr_cond)
        
            prev_block = find(p.block_order(1:n_block-1) == curr_cond, 1, 'last');
            prev_correct = behav_data.correct(:, prev_block);
            prev_lvl_order = p.trial_events(:, level_order_col, prev_block);

            probe_offsets = p.probe_offsets(curr_cond,:);
            step_size = squeeze(mean(p.staircases.step_size(:,end,:,curr_cond),1));

            new_probe_offsets = update_probe_offset(prev_correct, prev_lvl_order, probe_offsets, step_size, p.staircases.min_probe_offset, p.staircases.max_probe_offset);
            
            if p.disp_on
                disp('Old probe offsets:');
                disp(num2str(probe_offsets));
                disp('New probe offsets:');
                disp(num2str(new_probe_offsets)); 
            end

            p.probe_offsets(curr_cond,:) = new_probe_offsets;
        
            curr_probe_offsets = p.probe_offsets(curr_cond, p.trial_events(:, level_order_col, n_block))';
            probe_orientation = calc_probe_orientation(p.trial_events(:,test_orientation_col, n_block), curr_probe_offsets);
            p.trial_events(:, probe_orientation_col, n_block) = probe_orientation;
            
        end

    else

        probe_orientation = p.trial_events(:, probe_orientation_col, n_block);
    
    end

    corrected_probe_orientations = correct_orientation(probe_orientation)'; % Transform probe orientation to match the compass axis
    
    % Storing correct response
    cclockwise_trials = double(probe_orientation < p.trial_events(:,test_orientation_col, n_block));
    cclockwise_trials(cclockwise_trials == 0) = 2;
    p.correct_response(:,n_block) = cclockwise_trials;

    % Grab current block info
    n_trial = 1;

    % Task 
    trials_loop

    % Enter rest period
    if n_block < p.num_blocks

        % Briefly show progress
        rest_period_text = ['Block ' num2str(n_block) ' of ' num2str(p.num_blocks) 'completed'];
        rest_period_text_boundary = Screen('TextBounds', w.window, rest_period_text);
        rest_period_text_patch = CenterRectOnPoint(rest_period_text_boundary, w.centerX, w.centerY);
        Screen('DrawText', w.window, rest_period_text, rest_period_text_patch(1),  rest_period_text_patch(2), w.white);
        Screen('Flip', w.window);
        WaitSecs(2);

        % Enter rest period
        rest_period
    end

end

disp('Performance:')
for cond = 1:p.num_conds
    for lvl = 1:p.num_levels
        curr_trials = behav_data.correct(p.trial_events(:,3) == lvl, p.block_order == cond);
        behav_data.performance(lvl, cond) = mean(curr_trials(:));
        disp(['Condition ' p.cond_names{cond} ' Level ' num2str(lvl) ': ' num2str(round(100*behav_data.performance(lvl, cond))) '%']);
    end
end

%%

t.exp_end_time = GetSecs;
t.exp_dur = (t.exp_end_time - t.exp_start_time) / 60; % min
