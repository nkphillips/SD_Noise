%%% update_staircase

%% 2-down 1-up rule

% The current direction is stored as the previous direction
staircases.direction(curr_sc, 1, curr_sf) = staircases.direction(curr_sc, 2, curr_sf);

if staircases.responses(curr_sc, curr_sc_trial, curr_sf) == 1
    
    % Increase correct response counter
    staircases.correct_count(curr_sc, curr_sf) = staircases.correct_count(curr_sc, curr_sf) + 1;

    if p.demo_run
        disp(['Correct responses: ' num2str(staircases.correct_count(curr_sc, curr_sf))])
    end

    % Update contrast if there have been 2 correct responses
    if staircases.correct_count(curr_sc, curr_sf) == 2
        
        % Update direction
        staircases.direction(curr_sc, 2, curr_sf) = -1; % the current direction is updated as decreasing

        % Decrease contrast delta for next trial of the current staircase
        if curr_sc_trial < staircases.max_trials_per_sc
            staircases.contrast_deltas(curr_sc, curr_sc_trial+1, curr_sf) = max(staircases.min_contrast_delta, staircases.contrast_deltas(curr_sc, curr_sc_trial, curr_sf) - staircases.step_size(curr_sc, curr_sc_trial, curr_sf));
        end

        % Reset counter
        staircases.correct_count(curr_sc, curr_sf) = 0;
    
    elseif staircases.correct_count(curr_sc, curr_sf) == 1

        % Update direction
        staircases.direction(curr_sc, 2, curr_sf) = -1;
        
        % Next contrast delta is the same as the current contrast delta
        if curr_sc_trial < staircases.max_trials_per_sc
            staircases.contrast_deltas(curr_sc, curr_sc_trial+1, curr_sf) = staircases.contrast_deltas(curr_sc, curr_sc_trial, curr_sf);
        end
    end

else

    % Reset correct response counter
    staircases.correct_count(curr_sc, curr_sf) = 0;

    if p.demo_run
        disp('Correct response counter reset!')
    end

    % Update direction
    staircases.direction(curr_sc, 2, curr_sf) = 1; % the current direction is updated as increasing

    % Increase contrast delta for next trial of the current staircase
    if curr_sc_trial < staircases.max_trials_per_sc
        staircases.contrast_deltas(curr_sc, curr_sc_trial+1, curr_sf) = min(staircases.max_contrast_delta, staircases.contrast_deltas(curr_sc, curr_sc_trial, curr_sf) + staircases.step_size(curr_sc, curr_sc_trial, curr_sf));
    end

end

if p.demo_run && curr_sc_trial < staircases.max_trials_per_sc
    if staircases.correct_count(curr_sc, curr_sf) == 0
        disp(['Next contrast delta increased to ' num2str(100*staircases.contrast_deltas(curr_sc, curr_sc_trial+1, curr_sf)) '%'])
    elseif staircases.correct_count(curr_sc, curr_sf) == 2
        disp(['Next contrast delta decreased to ' num2str(100*staircases.contrast_deltas(curr_sc, curr_sc_trial+1, curr_sf)) '%'])
    end
end

%% Check for reversal

% If the direction of the current trial is different from the previous trial
if staircases.direction(curr_sc, 1, curr_sf) ~= staircases.direction(curr_sc, 2, curr_sf) && curr_sc_trial > 2

    % Increment reversal counter
    staircases.num_reversals(curr_sc, curr_sf) = staircases.num_reversals(curr_sc, curr_sf) + 1;

    % Record reversal index
    staircases.reversal_indices(curr_sc, curr_sc_trial, curr_sf) = 1;

    disp('Reversal detected.')

end

%% Update step size

% Reduce step size after some reversals
if staircases.num_reversals(curr_sc, curr_sf) >= staircases.max_reversals/2 && curr_sc_trial < staircases.max_trials_per_sc
    
    % Get the current step size
    old_step_size = staircases.step_size(curr_sc, curr_sc_trial, curr_sf);

    % Reduce step size
    staircases.step_size(curr_sc, curr_sc_trial+1, curr_sf) = max(staircases.min_step_size, staircases.step_size(curr_sc, curr_sc_trial, curr_sf) / 2);

    % Ensure step size does not exceed max step size
    staircases.step_size(curr_sc, curr_sc_trial+1, curr_sf) = min(staircases.max_step_size, staircases.step_size(curr_sc, curr_sc_trial+1, curr_sf));

    if p.demo_run
        disp(['Staircase step size reduced from ' num2str(100*old_step_size) ' to ' num2str(100*staircases.step_size(curr_sc, curr_sc_trial+1, curr_sf)) '%'])
    end

else

    % Next trial step size is the same as the current trial step size
    staircases.step_size(curr_sc, curr_sc_trial+1, curr_sf) = staircases.step_size(curr_sc, curr_sc_trial, curr_sf);

end

