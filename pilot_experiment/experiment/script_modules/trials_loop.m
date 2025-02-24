%%% trials_loop
% Main task loop.
% For every trial:
% Present test and reference SF (500-ms bandpass filtered noise; noise sample refreshed at 20 Hz),
% Record a response (2-s blank screen),
% Present ITI (blank screen) or top-up adaptation.
% Repeat missed trials.

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

redo_trial = 0;

for n_trial = 1:size(p.trial_events,1)

    % Get current test orientation
    curr_test_orient = p.trial_events(n_trial, 1, n_block); 

    % Create rotated probe line
    curr_probe_orient = p.trial_events(n_trial, 2, n_block); % Get current probe orientation
    curr_probe_orient_rad = deg2rad(curr_probe_orient);
    curr_rotation = [cos(curr_probe_orient_rad), -sin(curr_probe_orient_rad); ...
        sin(curr_probe_orient_rad), cos(curr_probe_orient_rad)];

    curr_probe_line = curr_rotation * (stimuli.probe_line_base - [w.centerX; w.centerY]);
    curr_probe_line = curr_probe_line + [w.centerX; w.centerY];

    % Get current contrast and filter width
    if curr_cond == 1

        curr_contrast = p.trial_events(n_trial, 3, n_block);
        curr_filter_width = 1;

    elseif curr_cond == 2

        curr_contrast = 1; 
        curr_filter_width = p.trial_events(n_trial, 3, n_block);

    end

    if p.demo_run
        disp(['Trial ' num2str(n_trial)])
        disp(['Test Orientation: ' num2str(curr_test_orient) '°'])
        disp(['Test Contrast: ' num2str(round(100*stimuli.contrast(curr_contrast),2)) '%'])
        disp(['Test Filter Width: ' num2str(stimuli.bp_filter_width(curr_filter_width)) '°'])
        if curr_probe_orient > 90
            disp(['Old Probe Orientation: ' num2str(curr_probe_orient-90) '°'])
        else
            disp(['Old Probe Orientation: ' num2str(curr_probe_orient+270) '°'])
        end
        disp(['Corrected Probe Orientation: ' num2str(curr_probe_orient) '°'])
    end

    %% Test orientation
    
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
            KbWait;
        end

    end
    
    %% Mask
    % present rapidly updating white noise
    
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
    
    %% Delay
    % blank period of a specified duration
    % just draw fixation

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
    
    %% Probe
    
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
            KbWait;
        end

    end

    %% Response Period

    % Initialize response state as "no response" being true
    no_response = 1;

    response_start = GetSecs;

    while no_response

        % Draw fixation
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

        % Flip
        Screen('Flip', w.window);

        % Check for response
        if ~p.simulate_response

            [key_pressed, first_press] = KbQueueCheck(p.device_number);
            which_press = find(first_press);

        else
            key_pressed = 1;
            which_press = datasample(p.keypress_numbers(1:2), 1); % random response
        end

        % If key is pressed
        if key_pressed && no_response

            % Find overlap between pressed keys and relevant keys
            relevant_keys = intersect(which_press, p.keypress_numbers(1:2));

            if ~isempty(relevant_keys) % If a relevant key was pressed
                
                first_relevant_key = relevant_keys(1); % Get the first relevant key

                % Store response based on the pressed key
                if first_relevant_key == p.keypress_numbers(1)
                  
                    behav_data.response(n_trial, n_block) = 1; % CCW
                    if p.demo_run
                        disp('Response: Test perceived as CCW.');
                    end
               
                elseif first_relevant_key == p.keypress_numbers(2)
                
                    behav_data.response(n_trial, n_block) = 2; % CW
                    if p.demo_run
                        disp('Response: Test perceived as CW.');
                    end
                
                end

                % Mark response as successful
                no_response = 0;

                if p.demo_run && p.correct_response(n_trial,n_block) == behav_data.response(n_trial, n_block)
                    disp('Correct response!')
                elseif p.demo_run && p.correct_response(n_trial,n_block) ~= behav_data.response(n_trial, n_block)
                    disp('Incorrect response!')
                end

            elseif any(which_press == p.keypress_numbers(end)) % Escape key
                
                if p.demo_run
                    disp('Escape key pressed. Exiting...');
                end
               
                return; % Exit the experiment
            
            end

        end

    end

    %% ITI

    if n_trial <= p.num_trials_per_block

        for n_frame = 1:frames.iti_frames_count(n_trial)

            % Draw fixation
            Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

            
            % Update frame deadlines
            if n_frame == 1
                iti_frames_onsets = frames.iti_frames_onsets{n_trial} + GetSecs;
            end

            % Flip
            [exe_timing.iti_VBLTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_StimulusOnsetTime{n_block}(n_trial,n_frame), exe_timing.iti_FlipTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_Missed{n_block}(n_trial,n_frame)] = ...
                Screen('Flip', w.window, iti_frames_onsets(n_frame));

        end

    end

    if p.demo_run, disp('  '); end

end