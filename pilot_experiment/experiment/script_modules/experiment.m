% Serial dependence + noise

%{


============================
Written by Luis D. Ramirez & Natalya Phillips
lur003@ucsd.edu
UCSD
Created 10.24.2024

%}

%%

if ~p.demo_run, wait_for_trigger; end
t.exp_start_time = GetSecs;

%% Block loop

for n_block = 1:p.num_blocks

    if p.demo_run, disp(['Block ' num2str(n_block)]); end

    % Pull the current block condition
    curr_cond = p.block_order(n_block);

    % Update probe offsets if it's not the first block of the condition
    if n_block > first_block(curr_cond)
    
        update_staircase
        
    else

        curr_probe_offsets = staircases.final_probe_offsets(curr_cond,:);

    end

    % Update probe orientations
    probe_orientation = p.trial_events(:,1,n_block) + curr_probe_offsets(p.trial_events(:,3,n_block)) .* datasample([-1 1], size(p.trial_events,1)); % test_orientation ± probe_offset
    p.trial_events(:,2,n_block) = correct_orientation(probe_orientation); % Transform probe orientation to match the compass axis

    % Grab current block info
    n_trial = 1;

    % Task 
    trials_loop

    % Enter rest period
    if n_block < p.num_blocks
        rest_period
    end

end

%%

t.exp_end_time = GetSecs;
t.exp_dur = (t.exp_end_time - t.exp_start_time) / 60; % min
