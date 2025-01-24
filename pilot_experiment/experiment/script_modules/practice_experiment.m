%%% practice_experiment


%%

if ~p.demo_run, wait_for_trigger; end
t.exp_start_time = GetSecs;

%% Block loop

for n_block = 1:p.num_blocks

    if p.demo_run, disp(['Block ' num2str(n_block)]); end

    % Grab current block info
    curr_adaptor = p.adapt_cond_order(1, n_block);
    curr_hemifield = p.adapt_cond_order(2, n_block);
    n_trial = 1;

    % Cue hemifield
    cue_hemifield

    % Task 
    practice_trials_loop

    % Enter rest period
    if n_block < p.num_blocks
        rest_period
    end

end

%%

t.exp_end_time = GetSecs;
t.exp_dur = (t.exp_end_time - t.exp_start_time) / 60; % min