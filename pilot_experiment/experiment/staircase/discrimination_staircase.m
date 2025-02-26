%%% detection_staircase
% 2 staircases per inducer SF that begin at the min and max contrast value

for curr_sf = 1:length(p.inducer_sfs)
   
    if p.demo_run
        disp(['SF = ' num2str(p.inducer_sfs(curr_sf)) ' cpd']);
    end

    for n_trial = 1:length(staircases.trial_order(:, curr_sf))
      
        %% Get current staircase info

        curr_sc = staircases.trial_order(n_trial, curr_sf);
        
        % Skip current trial if the max reversals have already been reached for the current staircase
        if staircases.num_reversals(curr_sc, curr_sf) == staircases.max_reversals
            % continue
        end
        
        % Get the current staircase trial number
        curr_sc_trial = find(n_trial == staircases.trial_indices(curr_sc, :, curr_sf));

        if p.demo_run   
            disp(['Staircase ' num2str(curr_sc) ', trial ' num2str(curr_sc_trial)])
        end

        %% Make contrast-adjusted textures

        % Get current contrast
        curr_contrast = min(1, stimuli.base_contrast + staircases.contrast_deltas(curr_sc, curr_sc_trial, curr_sf));

        % Display current contrast
        if p.demo_run
            disp(['Contrast delta = ' num2str(100*(curr_contrast - stimuli.base_contrast)) '%'])
        end

        % Store current contrast
        staircases.contrast(curr_sc, curr_sc_trial, curr_sf) = curr_contrast;

        % Make contrast-adjusted textures
        contrast_textures = squeeze(stimuli.full_contrast_textures(:,:, curr_sf, :)) .* curr_contrast;
        contrast_textures_made = nan(1, p.num_samples);
        for i = 1:p.num_samples
            contrast_textures_made(i) = Screen('MakeTexture', w.window, contrast_textures(:,:,i));
        end

        %% Draw textures

        n_noise_sample = 0;
        
        for n_frame = 1:frames.frames_count

            % Update noise sample
            if frames.noise_sample_update(n_frame)
                n_noise_sample = n_noise_sample + 1;
                curr_noise_sample = frames.update_noise_sample_seq(n_trial, n_noise_sample, curr_sf);
            end

            % Draw stimulus with correct contrast 
            if frames.contrast_delta(n_trial, n_frame, curr_sf) == 1
                Screen('DrawTexture', w.window, contrast_textures_made(curr_noise_sample))
            else
                Screen('DrawTexture', w.window, stimuli.textures_made(curr_sf, curr_noise_sample))
            end

            % Draw Aperture
            Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch)

            % Draw Fixation
            Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
            Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

            % Set frame onsets if the first frame
            if n_frame == 1
                frame_onsets = frames.frame_onsets + GetSecs;
            end

            % Flip stimuli
            Screen('Flip', w.window, frame_onsets(n_frame));

        end

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
                staircases.responses(curr_sc, curr_sc_trial, curr_sf) = 0;
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
                    staircases.responses(curr_sc, curr_sc_trial, curr_sf) = 1;
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
                    staircases.responses(curr_sc, curr_sc_trial, curr_sf) = 0;
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
                        staircases.responses(curr_sc, curr_sc_trial, curr_sf) = 1;
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

        if p.demo_run

            if no_response
                disp('Change missed')
            end

        end
        
        %% Update staircase

        update_staircase

        if p.demo_run
            disp(' ')
            % KbWait;
        end
        
    end
end