%% Loading complete screen

loading_text = 'Loading complete!';
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Start screen

if p.which_setup ~= 2
    
    start_screen_text = ['Press ' KbName(p.trigger_key) ' to begin.'];
    start_screen_text_boundary = Screen('TextBounds', w.window, start_screen_text);
    start_screen_text_patch = CenterRectOnPoint(start_screen_text_boundary, w.centerX, w.centerY);
    Screen('DrawText', w.window, start_screen_text, start_screen_text_patch(1),  start_screen_text_patch(2), w.white);
    
else
    
    DrawFormattedText(w.window, '~', w.ppd, w.ppd, w.white, w.bg_color);
    
end

Screen('Flip', w.window);

disp('Waiting for trigger...')

%% Wait for trigger

trigger_pressed = 0;

if ~p.simulate_response
    
    while ~trigger_pressed
        
        % Check for response
        [key_pressed, ~, first_press] = KbCheck(p.device_number);
        which_press = find(first_press);
        
        if key_pressed && which_press(1) == p.trigger_key
            trigger_pressed = 1;
        end
        
    end
    
    t.trigger_time_stamp = GetSecs;
    
end

disp('Trigger detected!')

%% Draw Fixation

Screen('DrawTexture', w.window, fixation_space_made, [], fixation_space_patch); % Fixation circle
Screen('FillOval', w.window, p.fixation_dot_color, fixation_dot_patch); % Fixation dot
Screen('Flip', w.window);
WaitSecs(2);