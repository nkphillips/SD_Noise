%%% detection_staircase
% 2 staircases per inducer SF that begin at the min and max contrast value

n_block = 0;

for curr_cond = 1:p.num_conds
    for curr_lvl = 1:p.num_levels

        n_block = n_block + 1;

        for n_trial = 1:length(staircases.trial_order(:, curr_lvl, curr_cond))
        
            %% Get current staircase info

            % Get current test orientation
            curr_test_orient = staircases.test_orientation(n_trial,curr_lvl,curr_cond);           
            
            % Get current staircase
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
            curr_probe_offset = staircases.probe_offsets(curr_sc,curr_sc_trial,curr_lvl,curr_cond);
            curr_probe_orientation = calc_probe_orientation(curr_test_orient,curr_probe_offset);

            % Create rotated probe line
            curr_probe_orientation_rad = deg2rad(curr_probe_orientation);
            curr_rotation = [cos(curr_probe_orientation_rad), -sin(curr_probe_orientation_rad); ...
                sin(curr_probe_orientation_rad), cos(curr_probe_orientation_rad)];

            curr_probe_line = curr_rotation * (stimuli.probe_line_base - [w.centerX; w.centerY]);
            curr_probe_line = curr_probe_line + [w.centerX; w.centerY];

            % Get current contrast and filter width
            if curr_cond == 1

                curr_contrast = curr_lvl;
                curr_filter_width = 1;

            elseif curr_cond == 2

                curr_contrast = 1;
                curr_filter_width = curr_lvl;

            end


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


            %% Draw Test

            n_noise_sample = 0;

            for n_frame = 1:frames.test_frames_count

                % Update noise sample when a new sample is needed
                if frames.test_noise_sample_update(n_frame)
                    n_noise_sample = n_noise_sample + 1;
                    curr_noise_sample = frames.test_noise_sample_update_seq(n_trial, n_noise_sample, n_block);
                end

                % Draw Test
                Screen('DrawTexture', w.window, stimuli.test_textures_made(curr_contrast, curr_filter_width, curr_noise_sample), [], noise_patch, curr_test_orient);

                % Stimulus aperture
                Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch, curr_test_orient);

                % Draw fixation
                Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
                Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

                % Flip
                if n_frame == 1
                    test_frames_onsets = frames.test_frames_onsets + GetSecs;
                end
                Screen('Flip', w.window, test_frames_onsets(n_frame));

                if p.demo_run && n_frame == frames.test_frames_count
                    % KbWait;
                end

            end

            %% Draw mask

            n_noise_sample = 0;

            for n_frame = 1:frames.mask_frames_count

                % Update noise sample when a new sample is needed
                if frames.mask_noise_sample_update(n_frame)
                    n_noise_sample = n_noise_sample + 1;
                    curr_noise_sample = frames.mask_noise_sample_update_seq(n_trial, n_noise_sample, n_block);
                end

                % Draw mask
                Screen('DrawTexture', w.window, stimuli.mask_textures_made(curr_contrast, curr_noise_sample), [], noise_patch);

                % Stimulus aperture
                Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch);

                % Draw fixation
                Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
                Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

                % Flip
                if n_frame == 1
                    mask_frames_onsets = frames.mask_frames_onsets + GetSecs;
                end
                Screen('Flip', w.window, mask_frames_onsets(n_frame)); % every frame has a deadline

                % if p.demo_run && n_frame == frames.test_frames_count
                %     KbWait;
                % end

            end

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

            %% Draw Probe

            for n_frame = 1:frames.probe_frames_count

                % Draw Line
                Screen('DrawLines', w.window, curr_probe_line, stimuli.probe_thickness, stimuli.probe_color);

                % Stimulus aperture
                Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch);

                % Draw fixation
                Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
                Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

                % Flip
                if n_frame == 1
                    probe_frames_onsets = frames.probe_frames_onsets + GetSecs;
                end
                Screen('Flip', w.window, probe_frames_onsets(n_frame));

                if p.demo_run && n_frame == frames.test_frames_count
                    %KbWait;
                end

            end

            %% Response

            % Initialize response state and response start
            no_response_recorded = 1;
            response_start = GetSecs;

            % Check for a key press
            if ~p.simulate_response

                [key_pressed, first_press] = KbQueueCheck(p.device_number);
                which_press = find(first_press);
                response_dur = GetSecs - response_start;

            else

                % Define psychometric function parameters
                threshold = 5; % Discrimination threshold in degrees
                slope = 0.3; % Controls steepness of psychometric function
                guess_rate = 0.5; % Chance performance (50% for 2AFC)
                lapse_rate = 0.02; % Small probability of lapses/mistakes
                
                % Calculate discrimination probability using cumulative normal function
                z_score = (curr_probe_offset - threshold) / slope;
                p_correct = guess_rate + (1 - guess_rate - lapse_rate) * normcdf(z_score);
                
                if p.demo_run
                    disp(['Discrimination probability: ' num2str(round(100*p_correct)) '%'])
                end
                
                % Simulate response based on discrimination probability
                key_pressed = rand() < p_correct;

                if key_pressed
                    which_press = p.keypress_numbers(1);
                end
            
                % Record response duration
                response_dur = GetSecs - response_start;

            end

            % Evaluate key press
            if key_pressed 

                % Find overlap between pressed keys and relevant keys
                relevant_keys = intersect(which_press, p.keypress_numbers(1));
                if ~isempty(relevant_keys) % If a relevant key was pressed

                    first_relevant_key = relevant_keys(1); % Get the first relevant key

                    % Store response based on the pressed key
                    if first_relevant_key == p.keypress_numbers(1)
                        staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 1;
                        staircases.response_dur(curr_sc, curr_sc_trial, curr_lvl) = response_dur;
                    end

                    % No response is now false
                    no_response_recorded = 0;

                elseif any(which_press == KbName('ESCAPE')) % Escape key

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
                    if ~p.simulate_response && no_response_recorded
                        [key_pressed, first_press] = KbQueueCheck(p.device_number);
                        which_press = find(first_press);
                        response_dur = GetSecs - response_start;
                    end

                    % Evaluate response
                    if key_pressed && no_response_recorded

                        % Find overlap between pressed keys and relevant keys
                        relevant_keys = intersect(which_press, p.keypress_numbers(1));

                        if ~isempty(relevant_keys) % If a relevant key was pressed

                            first_relevant_key = relevant_keys(1); % Get the first relevant key

                            % Store response based on the pressed key
                            if first_relevant_key == p.keypress_numbers(1)
                                staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 1;
                                staircases.response_dur(curr_sc, curr_sc_trial, curr_lvl) = response_dur;
                                if p.demo_run
                                    disp('Change detected')
                                end
                            end

                            % No response is now false
                            no_response_recorded = 0;

                        elseif any(which_press == KbName('ESCAPE')) % Escape key

                            if p.demo_run
                                disp('Escape key pressed. Exiting...');
                            end
                            return; % Exit the experiment

                        end

                    end

                end

                if no_response_recorded
                    staircases.responses(curr_sc, curr_sc_trial, curr_lvl) = 0;
                    if p.demo_run
                        disp('Change missed')
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
        
        rest_period

    end
end