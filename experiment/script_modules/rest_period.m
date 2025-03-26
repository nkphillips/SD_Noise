%%% rest_period
% Provide rest for a minimum amount of time
% Subject can proceed when ready

if p.simulate_response
    exit_rest = 1;
else
    exit_rest = 0;
end

%% Mandatory rest

if ~p.simulate_response
    for n_frame = 1:frames.rest_frames_count

        % Draw Fixation
        Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
        Screen('FillOval', w.window, p.rest_fixation_color, fixation_dot_patch); % Fixation dot

        % Flip
        % Update frame deadlines
        if n_frame == 1
            pres_timing.rest_flip_times{n_block} = frames.rest_frame_deadlines + GetSecs;
        end
        Screen('Flip', w.window, pres_timing.rest_flip_times{n_block}(n_frame));

    end
end

%% Optional rest

tic;
while ~exit_rest

    % Draw Fixation
    Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
    Screen('FillOval', w.window, w.blue, fixation_dot_patch); % Fixation dot

    % Flip
    Screen('Flip', w.window, pres_timing.rest_flip_times{n_block}(n_frame));

    % Check for response
    [key_pressed, first_press] = KbQueueCheck(p.device_number);
    which_press = find(first_press);

    % If key is pressed
    if key_pressed

        if which_press(1) == p.keypress_numbers(3)
            exit_rest = 1;
            t.optional_rest_dur(n_block) = toc;
        elseif which_press(1) == p.keypress_numbers(end)
            return;
        end
    end

end

% Draw fixation
Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot
Screen('Flip', w.window);

WaitSecs(2);
