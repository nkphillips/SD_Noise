function staircases = update_step_size(staircases, curr_sc, curr_lvl, curr_cond, curr_sc_trial, p)

% Reduce step size after some reversals
if staircases.num_reversals(curr_sc, curr_lvl, curr_cond) >= staircases.num_reversals_to_consider && curr_sc_trial < staircases.max_trials_per_sc
    
    % Get the current step size
    old_step_size = staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond);

    % Reduce step size
    staircases.step_size(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond) = max(staircases.min_step_size, staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond) / 2);

    if p.disp_on
        disp(['Staircase step size for next trial reduced from ' num2str(old_step_size) ' to ' num2str(staircases.step_size(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond))])
    end
    
end

end 