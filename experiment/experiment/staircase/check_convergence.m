function staircases = check_convergence(staircases, curr_sc, curr_lvl, curr_cond, p)

    % Find reversal indices
    reversal_idx = find(staircases.reversal_indices(curr_sc,:,curr_lvl,curr_cond) == 1);
    
    % Only check convergence if we have enough reversals
    if length(reversal_idx) >= staircases.num_reversals_to_consider
        
        % Get the last N reversals
        last_reversals_idx = reversal_idx(end-staircases.num_reversals_to_consider+1:end);
        
        % Get the trials between first and last of these reversals
        trial_range = last_reversals_idx(1):last_reversals_idx(end);
        
        % Calculate mean probe offset for these reversals
        mean_offset = mean(staircases.probe_offsets(curr_sc, last_reversals_idx, curr_lvl, curr_cond));
        
        % Find trials near this probe offset (within ±10% of mean probe offset)
        offset_window = 0.1 * mean_offset;
        offsets = staircases.probe_offsets(curr_sc, trial_range, curr_lvl, curr_cond);
        responses = staircases.responses(curr_sc, trial_range, curr_lvl, curr_cond);
        relevant_trials = abs(offsets - mean_offset) <= offset_window;
        
        % Calculate discrimination performance
        performance = mean(responses(relevant_trials)) * 100;
        
        % Check if performance is within window of expected 70.7%
        performance_window = 2; % percentage points
        target_performance = 70.7;
        
        if abs(performance - target_performance) <= performance_window
            staircases.converged(curr_sc, curr_lvl, curr_cond) = 1;
            
            if p.demo_run
                disp(['Staircase converged! Performance = ' num2str(performance, '%.1f') '%']);
            end
            
        end
    end
end 