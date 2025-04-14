%%% update_staircase

%% 2-down 1-up rule

% Get current probe offset and step size
curr_probe_offset = staircases.probe_offsets(curr_sc, curr_sc_trial, curr_lvl, curr_cond);
curr_step_size = staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond);

% The current direction is stored as the previous direction
staircases.direction(curr_sc, 1, curr_lvl, curr_cond) = staircases.direction(curr_sc, 2, curr_lvl, curr_cond);

% Initialize next trial's step size to current trial's step size
if curr_sc_trial < staircases.max_trials_per_sc
    staircases.step_size(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond) = staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond);
end


if staircases.responses(curr_sc, curr_sc_trial, curr_lvl, curr_cond) == 1
    
    % Increase correct response counter
    staircases.correct_count(curr_sc, curr_lvl, curr_cond) = staircases.correct_count(curr_sc, curr_lvl, curr_cond) + 1;

    if p.disp_on
        disp(['Correct responses: ' num2str(staircases.correct_count(curr_sc, curr_lvl, curr_cond))])
    end

    % If there have been 2 correct responses
    if staircases.correct_count(curr_sc, curr_lvl, curr_cond) == 2 

        % Update direction of current trial to decreasing
        staircases.direction(curr_sc, 2, curr_lvl, curr_cond) = -1; 

        % Check for reversal and possibly update step size
        if curr_sc_trial > 1 && staircases.direction(curr_sc, 1, curr_lvl, curr_cond) == 1

            % Increment reversal counter
            staircases.num_reversals(curr_sc, curr_lvl, curr_cond) = staircases.num_reversals(curr_sc, curr_lvl, curr_cond) + 1;

            % Record reversal index
            staircases.reversal_indices(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = 1;

            % Display message
            if p.disp_on
                disp('Reversal detected (increasing to decreasing).')
            end

            % Update step size first (affects next contrast delta)
            staircases = update_step_size(staircases, curr_sc, curr_lvl, curr_cond, curr_sc_trial, p);
            
            % Check convergence (might end staircase)
            staircases = check_convergence(staircases, curr_sc, curr_lvl, curr_cond, p);

        end

        % Now use the possibly updated step size to decrease contrast delta for next trial
        if curr_sc_trial < staircases.max_trials_per_sc

            % Display message
            if p.disp_on
                disp(['Current step size: ' num2str(staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond)) '°'])
            end

            % Decrease probe offset for next trial
            staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond) = round(max(staircases.min_probe_offset, staircases.probe_offsets(curr_sc, curr_sc_trial, curr_lvl, curr_cond) - staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond)));
            if p.disp_on
                disp(['Next probe offset decreased to ' num2str(staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond)) '°'])
            end

        end

        % Reset counter
        staircases.correct_count(curr_sc, curr_lvl, curr_cond) = 0;

    elseif staircases.correct_count(curr_sc, curr_lvl, curr_cond) == 1

        % Update direction
        staircases.direction(curr_sc, 2, curr_lvl, curr_cond) = -1;

        % Next probe offset is the same as the current probe offset
        if curr_sc_trial < staircases.max_trials_per_sc

            % Display message
            if p.disp_on
                disp(['Current step size: ' num2str(staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond)) '°'])
            end

            % Next probe offset is the same as the current probe offset
            staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond) = staircases.probe_offsets(curr_sc, curr_sc_trial, curr_lvl, curr_cond);
            if p.disp_on
                disp(['Next probe offset will stay as ' num2str(staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond)) '°'])
            end

        end

    end

else

    % Check for reversal and possibly update step size
    if curr_sc_trial > 1 && staircases.direction(curr_sc, 1, curr_lvl, curr_cond) == -1

        % Increment reversal counter
        staircases.num_reversals(curr_sc, curr_lvl, curr_cond) = staircases.num_reversals(curr_sc, curr_lvl, curr_cond) + 1;

        % Record reversal index
        staircases.reversal_indices(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = 1;

        % Display message
        if p.disp_on
            disp('Reversal detected (decreasing to increasing).')
        end

        % Update step size first (affects next contrast delta)
        staircases = update_step_size(staircases, curr_sc, curr_lvl, curr_cond, curr_sc_trial, p);
        
        % Check convergence (might end staircase)
        staircases = check_convergence(staircases, curr_sc, curr_lvl, curr_cond, p);

    end

    % Reset counter and update direction
    staircases.correct_count(curr_sc, curr_lvl, curr_cond) = 0;
    staircases.direction(curr_sc, 2, curr_lvl, curr_cond) = 1;

    % Now use the possibly updated step size to increase contrast delta for next trial
    if curr_sc_trial < staircases.max_trials_per_sc

        % Display message
        if p.disp_on
            disp(['Current step size: ' num2str(staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond)) '°'])
        end

        % Increase probe offset for next trial
        staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond) = round(min(staircases.max_probe_offset, staircases.probe_offsets(curr_sc, curr_sc_trial, curr_lvl, curr_cond) + staircases.step_size(curr_sc, curr_sc_trial, curr_lvl, curr_cond)));
        if p.disp_on
            disp(['Next probe offset increased to ' num2str(staircases.probe_offsets(curr_sc, curr_sc_trial+1, curr_lvl, curr_cond)) '°'])
        end
        
    end

end

if p.disp_on
    KbWait;
end
