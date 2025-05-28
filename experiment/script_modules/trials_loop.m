%%% trials_loop

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

for n_trial = 1:size(p.trial_events,1)

    curr_test_orient = p.trial_events(n_trial, test_orientation_col, n_block); 
    curr_probe_orient = p.trial_events(n_trial, probe_orientation_col, n_block);

    if curr_cond == 1
        curr_contrast = p.trial_events(n_trial, level_order_col, n_block);
        curr_filter_width = 1;

    elseif curr_cond == 2
        curr_contrast = 1; 
        curr_filter_width = p.trial_events(n_trial, level_order_col, n_block);

    end

    if p.disp_on
        disp(['Trial ' num2str(n_trial)])
        disp(['Test Contrast: ' num2str(round(100*p.contrast(curr_contrast),2)) '%'])
        disp(['Test Filter Width: ' num2str(p.orientation_bp_filter_width(curr_filter_width)) '°'])
        disp(['Test Orientation: ' num2str(round(curr_test_orient,2)) '°'])
        disp(['Probe Orientation: ' num2str(round(curr_probe_orient,2)) '°'])
    end
    
    %% Draw Test orientation
    
    n_noise_sample = 0;

    for n_frame = 1:frames.test_frames_count

        if frames.test_noise_sample_update(n_frame)
            n_noise_sample = n_noise_sample + 1;
            curr_noise_sample = frames.test_noise_sample_update_seq(n_trial, n_noise_sample, n_block);
        end

        Screen('DrawTexture', w.window, stimuli.test_textures_made(curr_contrast, curr_filter_width, curr_noise_sample), [], noise_patch, curr_test_orient);
        Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch, curr_test_orient);
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch);
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch);

        if n_frame == 1
           test_frames_onsets = frames.test_frames_onsets + GetSecs;
        end
        Screen('Flip', w.window, test_frames_onsets(n_frame));

        if p.disp_on && n_frame == frames.test_frames_count
            %pause(3);
        end

    end
    
    %% Draw Mask
    
    n_noise_sample = 0;

    for n_frame = 1:frames.mask_frames_count

        if frames.mask_noise_sample_update(n_frame)
            n_noise_sample = n_noise_sample + 1;
            curr_noise_sample = frames.mask_noise_sample_update_seq(n_trial, n_noise_sample, n_block);
        end

        Screen('DrawTexture', w.window, stimuli.mask_textures_made(curr_contrast, curr_noise_sample), [], noise_patch);
        Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch);
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); 
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); 

        if n_frame == 1
            mask_frames_onsets = frames.mask_frames_onsets + GetSecs;
        end
        Screen('Flip', w.window, mask_frames_onsets(n_frame));

        if p.disp_on && n_frame == frames.mask_frames_count
            %pause(3);
        end

    end
    
    %% Delay period

    for n_frame = 1:frames.delay_frames_count

        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); 
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); 

        if n_frame == 1
            delay_frames_onsets = frames.delay_frames_onsets + GetSecs;
        end
        Screen('Flip', w.window, delay_frames_onsets(n_frame));

        if p.disp_on && n_frame == frames.delay_frames_count
            %pause(3);
        end

    end
    
    %% Draw Probe

    KbQueueFlush(p.device_number,1);

    for n_frame = 1:frames.probe_frames_count

        Screen('DrawTexture',w.window, stimuli.probe_line_made, [], probe_patch, curr_probe_orient);
        Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch);
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch);
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch);

        if n_frame == 1
            probe_frames_onsets = frames.probe_frames_onsets + GetSecs;
        end
        Screen('Flip', w.window, probe_frames_onsets(n_frame));
        
        if p.disp_on && n_frame == frames.probe_frames_count
            %pause(3);
        end

    end

    %% Response Period

    no_response_recorded = 1;

    response_start = GetSecs;

    while no_response_recorded

        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch);
        Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); 

        Screen('Flip', w.window);

        if ~p.simulate_response
            [key_pressed, first_press] = KbQueueCheck(p.device_number);
            which_press = find(first_press);
            response_dur = GetSecs - response_start;
        else
            key_pressed = 1;
            which_press = datasample(p.keypress_numbers, 1);
            response_dur = 0;
        end

        if key_pressed && no_response_recorded

            relevant_keys = intersect(which_press, p.keypress_numbers);

            if ~isempty(relevant_keys)
                
                first_relevant_key = relevant_keys(1); 

                if first_relevant_key == p.keypress_numbers(1)
                    behav_data.response(n_trial, n_block) = 1; 
                    if p.disp_on, disp('Response: CCW'); end
                elseif first_relevant_key == p.keypress_numbers(2)
                    behav_data.response(n_trial, n_block) = 2;
                    if p.disp_on, disp('Response: CW'); end
                end

                no_response_recorded = 0;

                behav_data.correct(n_trial, n_block) = p.correct_response(n_trial, n_block) == behav_data.response(n_trial, n_block);
                behav_data.response_dur(n_trial, n_block) = response_dur;
                    
                if p.disp_on && behav_data.correct(n_trial, n_block)
                    disp('Correct response!')
                elseif p.disp_on && ~behav_data.correct(n_trial, n_block)
                    disp('Incorrect response!')
                end

            elseif any(which_press == KbName('ESCAPE')) 
                if p.disp_on, disp('Escape key pressed. Exiting...'); end
                break;
            
            end
       
        end

    end

    %% ITI

    if n_trial <= p.num_trials_per_block

        for n_frame = 1:frames.iti_frames_count(n_trial)

            Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch);

            if n_frame == 1
                iti_frames_onsets = frames.iti_frames_onsets{n_trial} + GetSecs;
            end

            [exe_timing.iti_VBLTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_StimulusOnsetTime{n_block}(n_trial,n_frame), exe_timing.iti_FlipTimestamp{n_block}(n_trial,n_frame), exe_timing.iti_Missed{n_block}(n_trial,n_frame)] = ...
                Screen('Flip', w.window, iti_frames_onsets(n_frame));

        end

    end

    if p.disp_on, disp('  '); end

end