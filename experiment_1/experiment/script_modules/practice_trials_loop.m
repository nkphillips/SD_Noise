%%% practice_trials_loop
% Main task practice loop.
% For every trial:
% Present test and reference SF (500-ms bandpass filtered noise; noise sample refreshed at 20 Hz),
% Record a response (2-s blank screen),
% Present ITI (blank screen)
% Repeat missed trials.


%{ 

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

redo_trial = 0;

while n_trial <= p.num_trials_per_block
    
    KbEventFlush(p.device_number);

    % Grab current trial info
    curr_test_sf_indx = find(p.test_sf_order(n_block,n_trial) == p.test_sfs);

    if p.demo_run
        disp(['Trial #' num2str(n_trial)])
        disp(['Test SF: ' num2str(round(p.sf_ratios(curr_test_sf_indx),2)) ' oct']);
    end

    % Initialize noise samples
    curr_test_noise_sample = randperm(p.num_test_samples, 1);
    curr_reference_noise_sample = randperm(p.num_reference_samples, 1);

    %% Present test and reference SF

    n_test_sample = 0;
    n_ref_sample = 0;

    for n_frame = 1:frames.test_frames_count

        % Update frame deadlines
        if n_frame == 1

            pres_timing.test_flip_times{n_block}(n_trial,:) = frames.test_frame_deadlines + GetSecs;

        end

        % Update noise sample
        if frames.test_noise_sample_update(n_frame)

            n_test_sample = n_test_sample + 1;
            n_ref_sample = n_ref_sample + 1;

            curr_test_noise_sample = frames.test_noise_sample_update_seq{n_block}(n_trial,n_test_sample);
            curr_reference_noise_sample = frames.reference_noise_sample_update_seq{n_block}(n_trial,n_ref_sample);

        end

        % Draw test and reference SF
        if curr_hemifield == find(strcmp(p.hemifields,'L'))

            Screen('DrawTexture', w.window, stimuli.test_sf_made(curr_test_sf_indx, curr_test_noise_sample), [], noise_L_patch)
            Screen('DrawTexture', w.window, stimuli.reference_sf_made(curr_reference_noise_sample), [], noise_R_patch)

        elseif curr_hemifield == find(strcmp(p.hemifields,'R'))

            Screen('DrawTexture', w.window, stimuli.reference_sf_made(curr_reference_noise_sample), [], noise_L_patch)
            Screen('DrawTexture', w.window, stimuli.test_sf_made(curr_test_sf_indx, curr_test_noise_sample), [], noise_R_patch)

        end

        % Draw Apertures
        Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch)

        % Draw Fixation
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

        % Flip
        [exe_timing.test_VBLTimestamp{n_block}(n_trial,n_frame), exe_timing.test_StimulusOnsetTime{n_block}(n_trial,n_frame), exe_timing.test_FlipTimestamp{n_block}(n_trial,n_frame), exe_timing.test_Missed{n_block}(n_trial,n_frame)] = ...
            Screen('Flip', w.window, pres_timing.test_flip_times{n_block}(n_trial,n_frame));

        % WaitSecs(0.250);

    end

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

            if which_press(1) == p.keypress_numbers(1) || which_press(1) == p.keypress_numbers(2) % If it's a task relevant key

                if p.demo_run
                    disp(' ')
                    disp([KbName(which_press(1)) ' response detected!'])
                end

                no_response = 0; % A relevant key has been pressed, so flip "no response" to false

                % Store response as whether the test is perceived as lower (0) / greater (1) SF than the reference
                if which_press(1) == p.keypress_numbers(1)

                    behav_data.response{n_block}(n_trial) = 0;

                    % Update fixation dot color for feedback
                    if p.test_sf_order(n_block,n_trial) <= p.reference_sf
                        p.fixation_dot_color = w.green;
                    else
                        p.fixation_dot_color = w.red;
                    end

                elseif which_press(1) == p.keypress_numbers(2)

                    behav_data.response{n_block}(n_trial) = 1;

                    % Update fixation dot color for feedback
                    if p.test_sf_order(n_block,n_trial) >= p.reference_sf
                        p.fixation_dot_color = w.green;
                    else
                        p.fixation_dot_color = w.red;
                    end

                end

            elseif  which_press(1) == p.keypress_numbers(end) % If escape is pressed, quit the experiment

                return;

            end
            
        end

    end

    %% Check if there was no response and setup appropriate vars to redo the missed trial
    % If there's no response, redo the previous and current trial
    % If the previous test SF is missed again, stay there
    % only proceed once both the previous and current trial get a response

    if no_response

        % Update trackers and counters
        behav_data.missed_trial_indices{n_block} = [behav_data.missed_trial_indices{n_block} n_trial];
        behav_data.num_trials_per_block(n_block) = behav_data.num_trials_per_block(n_block) + 1;

        if n_trial == 1

            behav_data.test_sf_order{n_block} = [behav_data.test_sf_order{n_block}(1), behav_data.test_sf_order{n_block}];

            if p.demo_run, disp('No response! Previous trial will repeat.'); end

        elseif redo_trial == 0

            behav_data.test_sf_order{n_block} = [behav_data.test_sf_order{n_block}(1:n_trial), behav_data.test_sf_order{n_block}(n_trial-1:end)];

            n_trial = n_trial - 1; % when going back a trial, stay there until it is not missed, do not go back even more. (note that this n-1 trial is irrelavant to the analysis)

            redo_trial = 1;

            if p.demo_run, disp('No response! Previous trial will repeat.'); end

        else

            behav_data.test_sf_order{n_block} = [behav_data.test_sf_order{n_block}(1:n_trial), behav_data.test_sf_order{n_block}(n_trial:end)];
            if p.demo_run, disp('No response! Redo trial will be reattempted.'); end

        end

    else

        if p.demo_run && redo_trial, disp('Redo successful.'); end

        n_trial = n_trial + 1;
        redo_trial = 0;

    end

    %% ITI

    n_sample = 0;

    if n_trial <= p.num_trials_per_block

        % Get initial noise sample
        curr_noise_sample = randperm(p.num_adaptor_samples, 1);

        for n_frame = 1:frames.top_up_adaptation_frames_count
            
            % Check for noise sample update
            if frames.top_up_adaptation_noise_sample_update(n_frame)

                n_sample = n_sample + 1;
                curr_noise_sample = frames.top_up_adaptation_noise_sample_update_seq{n_block}(n_trial,n_sample);

            end

            % Draw adaptor
            if curr_adaptor > 1
                if curr_hemifield == find(strcmp(p.hemifields,'L'))

                    Screen('DrawTexture', w.window, stimuli.adaptor_sf_made(curr_adaptor-1,curr_noise_sample), [], noise_L_patch);

                elseif curr_hemifield == find(strcmp(p.hemifields,'R'))

                    Screen('DrawTexture', w.window, stimuli.adaptor_sf_made(curr_adaptor-1,curr_noise_sample), [], noise_R_patch);

                end
            end

            % Draw aperture
            Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch)

            % Draw fixation
            Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
            Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot

            % Flip
            % Update frame deadlines
            if n_frame == 1
                pres_timing.top_up_flip_times{n_block}(n_trial,:) = frames.top_up_adaptation_frame_deadlines + GetSecs;
            end
            [exe_timing.top_up_VBLTimestamp{n_block}(n_trial,n_frame), exe_timing.top_up_StimulusOnsetTime{n_block}(n_trial,n_frame), exe_timing.top_up_FlipTimestamp{n_block}(n_trial,n_frame), exe_timing.top_up_Missed{n_block}(n_trial,n_frame)] = ...
                Screen('Flip', w.window, pres_timing.top_up_flip_times{n_block}(n_trial,n_frame));

        end

        % Draw fixation
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot
        Screen('Flip', w.window);
        WaitSecs(t.gap_after_adapt);
   
    end

    if p.demo_run, disp('  '); end

    p.fixation_dot_color = w.black;

end