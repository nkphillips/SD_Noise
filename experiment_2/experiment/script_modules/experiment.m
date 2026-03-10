% Serial dependence + noise

%{


============================
Written by Luis D. Ramirez & Natalya Phillips
lur003@ucsd.edu
UCSD
Created 10.24.2024

%}

%%

if ~p.simulate_response, wait_for_trigger; end
t.exp_start_time = GetSecs;

%% Block loop

for n_block = 1:p.num_blocks

    if p.disp_on, disp(['Block ' num2str(n_block)]); end

    % Pull the current block condition
    curr_feature = p.feature_order(n_block);

    probe_orientation = p.trial_events(:, probe_orientation_col, n_block);

    % Storing correct response
    cclockwise_trials = double(probe_orientation < p.trial_events(:,test_orientation_col, n_block));
    cclockwise_trials(cclockwise_trials == 0) = 2;
    p.correct_response(:,n_block) = cclockwise_trials;

    % Grab current block info
    n_trial = 1;

    % Task
    trials_loop

    % Enter rest period
    if n_block < p.num_blocks

        % Briefly show progress
        rest_period_text = ['Block ' num2str(n_block) ' of ' num2str(p.num_blocks) ' completed'];
        rest_period_text_boundary = Screen('TextBounds', w.window, rest_period_text);
        rest_period_text_patch = CenterRectOnPoint(rest_period_text_boundary, w.centerX, w.centerY);
        Screen('DrawText', w.window, rest_period_text, rest_period_text_patch(1),  rest_period_text_patch(2), w.white);
        Screen('Flip', w.window);
        WaitSecs(2);

        % Enter rest period
        rest_period
    end

end

disp('')
for feature = 1:p.num_features
    for lvl = 1:p.num_levels
        curr_trials = behav_data.correct(p.trial_events(:,3) == lvl, p.feature_order == feature);
        behav_data.performance(lvl, feature) = mean(curr_trials(:));
    end
end
disp(['Performance: ' num2str(round(100*mean(behav_data.performance(:)))) '%']);

t.exp_end_time = GetSecs;
t.exp_dur = (t.exp_end_time - t.exp_start_time) / 60; % min

%% Save textures

if save_textures

    screen_text = 'Saving textures...';
    screen_text_boundary = Screen('TextBounds', w.window, screen_text);
    screen_text_patch = CenterRectOnPoint(screen_text_boundary, w.centerX, w.centerY);

    Screen('DrawText', w.window, screen_text, screen_text_patch(1),  screen_text_patch(2), w.white);
    Screen('Flip', w.window);

    % Remove fields from stimuli struct associated with the output of Screen('MakeTexture',...)
    struct_fieldnames = fieldnames(stimuli);
    struct_fieldnames(contains(struct_fieldnames,'made')) = [];

    % Store stimulus into textures struct
    for n_field = 1:numel(struct_fieldnames)
        textures.(struct_fieldnames{n_field}) = stimuli.(struct_fieldnames{n_field});
    end

    % Enter directory for saving the textures
    cd(dirs.texture_dir)

    if ~exist(textures_filename, "file")
        save(textures_filename,'textures','-mat','-v7.3');
    end

    screen_text = 'Textures saved!';
    screen_text_boundary = Screen('TextBounds', w.window, screen_text);
    screen_text_patch = CenterRectOnPoint(screen_text_boundary, w.centerX, w.centerY);

    Screen('DrawText', w.window, screen_text, screen_text_patch(1),  screen_text_patch(2), w.white);
    Screen('Flip', w.window);

    cd(dirs.script_dir)

end