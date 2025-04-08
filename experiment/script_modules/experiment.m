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

    % Update probe offsets if it's not the first block of the current condition
    if n_block > first_block(curr_cond)
    
        %update_staircase
        
    else

        curr_probe_offsets = p.probe_offsets(p.trial_events(:, level_order_col, n_block), curr_cond);

    end

    % Calculate probe orientations (test_orientation ± probe_offset)
    probe_orientation = calc_probe_orientation(p.trial_events(:,test_orientation_col, n_block), curr_probe_offsets);
    p.trial_events(:,probe_orientation_col,n_block) = correct_orientation(probe_orientation); % Transform probe orientation to match the compass axis

    % Grab current block info
    n_trial = 1;

    % Task 
    trials_loop

    % Enter rest period
    if n_block < p.num_blocks

        % Briefly show progress
        rest_period_text = ['Block ' num2str(n_block) ' of ' num2str(p.num_blocks)];
        rest_period_text_boundary = Screen('TextBounds', w.window, rest_period_text);
        rest_period_text_patch = CenterRectOnPoint(rest_period_text_boundary, w.centerX, w.centerY);
        Screen('DrawText', w.window, rest_period_text, rest_period_text_patch(1),  rest_period_text_patch(2), w.white);
        Screen('Flip', w.window);
        WaitSecs(2);

        % Enter rest period
        rest_period
    end

end

behav_data.performance = mean(behav_data.correct(:));

%%

t.exp_end_time = GetSecs;
t.exp_dur = (t.exp_end_time - t.exp_start_time) / 60; % min
