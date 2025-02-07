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

while n_trial <= p.num_trials_per_block

    % Grab current trial info
    curr_test_orient = p.trial_events(n_trial, 1, n_block);
    curr_probe_orient = p.trial_events(n_trial, 2, n_block);
    
    if p.demo_run
        disp(['Trial #' num2str(n_trial)])
        disp(['Test Orientation: ' num2str(round(curr_test_orient,2)) '°']);
        disp(['Probe Orientation: ' num2str(round(curr_probe_orient,2)) '°']);
    end

    %% Test orientation
    
    for n_frame = 1:frames.target_frames_count
        
        
        % The following variables and functions must be used to display the test stimulus
        %{

        Variables:
           
            frames.test_noise_sample_update
            frames.test_frame_deadlines
            stimuli.

stimuli.textures_made
stimuli.aperture_made
fixation_space_made
            
        Functions:
            
            Screen('DrawTexture', w.window, []);
            Screen('Flip', w.window, []);

        %}



        
    end
    
    %% Mask
    % present rapidly updating white noise
    
    
    %% Delay
    % blank period of a specified duration
    % just draw fixation

    
    %% Probe
    
    

    %% Response Period

    % Initialize response state as "no response" being true
    no_response = 1;

    response_start = GetSecs;

    while GetSecs - response_start < t.response_dur

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
                  
                    behav_data.response{n_block}(n_trial) = 0; % CCW
                    if p.demo_run
                        disp('Response: Test perceived as CCW.');
                    end
               
                elseif first_relevant_key == p.keypress_numbers(2)
                
                    behav_data.response{n_block}(n_trial) = 1; % CW
                    if p.demo_run
                        disp('Response: Test perceived as CW.');
                    end
                
                end

                % Mark response as successful
                no_response = 0;

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

        for n_frame = 1:frames.iti_frames_count

            % Draw fixation
            Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
            Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

            % Flip
            % Update frame deadlines
            [exe_timing.iti_VBLTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_StimulusOnsetTime{n_block}(n_trial,n_frame), exe_timing.iti_FlipTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_Missed{n_block}(n_trial,n_frame)] = ...
                Screen('Flip', w.window, pres_timing.iti_flip_times{n_block}(n_trial,n_frame));

        end

        % Draw fixation
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot
        Screen('Flip', w.window);

    end

    if p.demo_run, disp('  '); end

end