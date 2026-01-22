%%% detection_staircase
% 2 staircases per condition

n_block = 0;

for curr_cond = 1:p.num_conds
    for curr_lvl = 1:p.num_levels

        n_block = n_block + 1;

        for n_trial = 1:length(staircases.trial_order(:, curr_lvl, curr_cond))

            KbQueueFlush(p.device_number);
            if p.disp_on, disp(' '); end

            %% Get current staircase info
            
            curr_sc = staircases.trial_order(n_trial, curr_lvl, curr_cond);
            
            if staircases.num_reversals(curr_sc, curr_lvl, curr_cond) == staircases.max_reversals
                % continue
            end

            curr_sc_trial = find(n_trial == staircases.trial_indices(curr_sc, :, curr_lvl, curr_cond));
            curr_test_orient = staircases.test_orientation(curr_sc_trial, curr_lvl, curr_cond);           
            curr_probe_offset = staircases.probe_offsets(curr_sc, curr_sc_trial, curr_lvl, curr_cond);
            curr_probe_orientation = calc_probe_orientation(curr_test_orient, curr_probe_offset);
        
            orient_diff = abs(curr_probe_orientation - curr_test_orient);
            is_CW = curr_probe_orientation > curr_test_orient;
            
            if p.disp_on
                
                disp(['Staircase ' num2str(curr_sc) ', trial ' num2str(curr_sc_trial)])
                disp(['Test Orientation: ' num2str(curr_test_orient) '°'])
                disp(['Probe Orientation: ' num2str(curr_probe_orientation) '°'])
                disp(['Orientation difference: ' num2str(orient_diff) '°'])
                disp(' ')
                
                if is_CW
                    disp('Correct response: Right')
                else
                    disp('Correct response: Left')
                end

                % disp(['Test Contrast: ' num2str(round(100*stimuli.contrast(curr_contrast),2)) '%'])
                % disp(['Test Filter Width: ' num2str(stimuli.bp_filter_width(curr_filter_width)) '°'])

            end

            curr_probe_orientation = correct_orientation(curr_probe_orientation);

            % Create rotated probe line
            curr_probe_orientation_rad = deg2rad(curr_probe_orientation);
            curr_rotation = [cos(curr_probe_orientation_rad), -sin(curr_probe_orientation_rad); ...
                sin(curr_probe_orientation_rad), cos(curr_probe_orientation_rad)];

            curr_probe_line = curr_rotation * (stimuli.probe_line_base - [w.centerX; w.centerY]);
            curr_probe_line = curr_probe_line + [w.centerX; w.centerY];

            if curr_cond == 1
                curr_contrast = curr_lvl;
                curr_filter_width = 1;

            elseif curr_cond == 2
                curr_contrast = 1;
                curr_filter_width = curr_lvl;

            end

            %% Draw Test

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
                    % KbWait;
                end

            end

            %% Draw mask

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

                if p.disp_on && n_frame == frames.test_frames_count
                    % KbWait;
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

            end

            %% Draw Probe

            for n_frame = 1:frames.probe_frames_count

                Screen('DrawLines', w.window, curr_probe_line, p.probe_thickness, p.probe_color);
                Screen('DrawTexture', w.window, stimuli.aperture_made, [], aperture_patch);
                Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch);
                Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch);

                if n_frame == 1
                    probe_frames_onsets = frames.probe_frames_onsets + GetSecs;
                end
                Screen('Flip', w.window, probe_frames_onsets(n_frame));

                if p.disp_on && n_frame == frames.probe_frames_count
                    % KbWait;
                end

            end

            %% Response

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
                    mu = 0; 
                    sigma = 1;
                    guess_rate = 0; 
                    
                    p_CW = calc_pCW(orient_diff, mu, sigma, guess_rate);
                    response = rand() < p_CW;
                    
                    if response
                        which_press = p.keypress_numbers(2);
                    else
                        which_press = p.keypress_numbers(1);
                    end
                    
                    key_pressed = true;
                    response_dur = GetSecs - response_start;

                end

                if key_pressed && no_response_recorded

                    relevant_keys = intersect(which_press, p.keypress_numbers);

                    if ~isempty(relevant_keys) 

                        first_relevant_key = relevant_keys(1); 

                        if first_relevant_key == p.keypress_numbers(1) && ~is_CW
                            staircases.responses(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = 1;
                            staircases.response_dur(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = response_dur;
                            if p.disp_on, disp('Response: CCW'); end
                        
                        elseif first_relevant_key == p.keypress_numbers(2) && is_CW
                            staircases.responses(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = 1;
                            staircases.response_dur(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = response_dur;
                            if p.disp_on, disp('Response: CW'); end
                        
                        else
                            staircases.responses(curr_sc, curr_sc_trial, curr_lvl, curr_cond) = 0;
                            if p.disp_on && is_CW
                                disp('Response: CCW'); 
                            elseif p.disp_on && ~is_CW
                                disp('Response: CW'); 
                            end
                        end

                        no_response_recorded = 0;

                    elseif any(which_press == KbName('ESCAPE')) 
                        if p.disp_on
                            disp('Escape key pressed. Exiting...');
                        end
                        return; 

                    end

                end
            end

            %% ITI

            if n_trial < length(staircases.trial_order(:, curr_lvl, curr_cond))

                for n_frame = 1:frames.iti_frames_count(n_trial)

                    Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch);
        
                    if n_frame == 1
                        iti_frames_onsets = frames.iti_frames_onsets{n_trial} + GetSecs;
                    end
        
                    Screen('Flip', w.window, iti_frames_onsets(n_frame));
        
                end

            end

            %% Update staircase

            update_staircase

        end


        %% Rest period
        
        rest_period

    end
end