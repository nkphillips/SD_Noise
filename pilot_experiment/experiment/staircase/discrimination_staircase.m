%%% detection_staircase
% 2 staircases per inducer SF that begin at the min and max contrast value

for curr_cond = 1:p.num_conds
    for curr_lvl = 1:p.num_levels

        for n_trial = 1:length(staircases.trial_order(:, curr_lvl, curr_cond))

            % Get current test orientation
            curr_test_orient = [];

            if p.demo_run
                disp(['Trial ' num2str(n_trial)])
                disp(['Test Orientation: ' num2str(curr_test_orient) '°'])
                % disp(['Test Contrast: ' num2str(round(100*stimuli.contrast(curr_contrast),2)) '%'])
                % disp(['Test Filter Width: ' num2str(stimuli.bp_filter_width(curr_filter_width)) '°'])
                % if curr_probe_orient > 90
                %     disp(['Old Probe Orientation: ' num2str(curr_probe_orient-90) '°'])
                % else
                %     disp(['Old Probe Orientation: ' num2str(curr_probe_orient+270) '°'])
                % end
                % disp(['Corrected Probe Orientation: ' num2str(curr_probe_orient) '°'])
            end

            %% Get current staircase info

            curr_sc = staircases.trial_order(n_trial, curr_lvl, curr_cond);

            % Skip current trial if the max reversals have already been reached for the current staircase
            if staircases.num_reversals(curr_sc, curr_lvl, curr_cond) == staircases.max_reversals
                % continue
            end

            % Get the current staircase trial number
            curr_sc_trial = find(n_trial == staircases.trial_indices(curr_sc, :, curr_lvl, curr_cond));

            if p.demo_run
                disp(['Staircase ' num2str(curr_sc) ', trial ' num2str(curr_sc_trial)])
            end


            % Probe offset
            curr_probe_orientation = [];

            % Create rotated probe line
            curr_probe_orientation_rad = deg2rad(curr_probe_orientation);
            curr_rotation = [cos(curr_probe_orientation_rad), -sin(curr_probe_orientation_rad); ...
                sin(curr_probe_orientation_rad), cos(curr_probe_orientation_rad)];

            curr_probe_line = curr_rotation * (stimuli.probe_line_base - [w.centerX; w.centerY]);
            curr_probe_line = curr_probe_line + [w.centerX; w.centerY];

            
            %% Draw Test


            %% Draw mask


            %% Delay period

            for n_frame = 1:frames.delay_frames_count

                % Draw fixation
                Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
                Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

                % Flip
                if n_frame == 1
                    delay_frames_onsets = frames.delay_frames_onsets + GetSecs;
                end
                Screen('Flip', w.window, delay_frames_onsets(n_frame)); % every frame has a deadline

            end

            %% Draw probe


            %% Response

            % Initialize response state and response start
            no_response = 1;
            response_start = GetSecs;

            % Check for a key press
            if ~p.simulate_response

                [key_pressed, first_press] = KbQueueCheck(p.device_number);
                which_press = find(first_press);
                response_dur = GetSecs - response_start;

            else

                key_pressed = datasample([0 1], 1, 'Weights',[1-presumed_target presumed_target]);
                if key_pressed
                    which_press = p.keypress_numbers(1);
                    response_dur = GetSecs - response_start;
                    if p.demo_run
                        disp('Change detected')
                    end
                else
                    staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 0;
                    if p.demo_run
                        disp('Change missed')
                    end
                end

            end

            % Evaluate response
            if key_pressed && no_response

                % Find overlap between pressed keys and relevant keys
                relevant_keys = intersect(which_press, [0 p.keypress_numbers(1)]);
                if ~isempty(relevant_keys) % If a relevant key was pressed

                    response_dur(n_trial) = response_dur; % Record response duration
                    first_relevant_key = relevant_keys(1); % Get the first relevant key

                    % Store response based on the pressed key
                    if first_relevant_key == p.keypress_numbers(1)
                        staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 1;
                    end

                    % No response is now false
                    no_response = 0;

                elseif any(which_press == p.keypress_numbers(end)) % Escape key

                    if p.demo_run
                        disp('Escape key pressed. Exiting...');
                    end
                    return; % Exit the experiment

                end

            end

            %% ITI

            if n_trial < length(staircases.trial_order(:, curr_lvl, curr_cond))
                iti_start = GetSecs;

                while GetSecs - iti_start <= t.iti_dur

                    % Draw Fixation
                    Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

                    % Flip stimuli
                    Screen('Flip', w.window);

                    % Check for a key press
                    if ~p.simulate_response
                        [key_pressed, first_press] = KbQueueCheck(p.device_number);
                        which_press = find(first_press);
                        response_dur = GetSecs - response_start;
                    else

                        key_pressed = datasample([0 1], 1, 'Weights',[1-presumed_target presumed_target]);
                        if key_pressed
                            which_press = p.keypress_numbers(1);
                            response_dur = GetSecs - response_start;
                        else
                            staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 0;
                        end

                    end

                    % Evaluate response
                    if key_pressed && no_response

                        % Find overlap between pressed keys and relevant keys
                        relevant_keys = intersect(which_press, [0 p.keypress_numbers(1)]);

                        if ~isempty(relevant_keys) % If a relevant key was pressed

                            response_dur(n_trial) = response_dur; % Record response duration
                            first_relevant_key = relevant_keys(1); % Get the first relevant key

                            % Store response based on the pressed key
                            if first_relevant_key == p.keypress_numbers(1)
                                staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 1;
                                if p.demo_run
                                    disp('Change detected')
                                end
                            end

                            % No response is now false
                            no_response = 0;

                        elseif any(which_press == p.keypress_numbers(end)) % Escape key

                            if p.demo_run
                                disp('Escape key pressed. Exiting...');
                            end
                            return; % Exit the experiment

                        end

                    end

                end
            end

            %% Update staircase

            update_staircase

            if p.demo_run
                disp(' ')
                % KbWait;
            end

        end


        %% Rest period

        

    end
end