%%% run_experiment

%{ 

Written by Luis D. Ramirez & Natalya Phillips
UCSD
lur003@ucsd.edu

%}


function run_info = run_experiment(p,w,dirs)

%% Verify run number to make sure nothing gets overwritten
% Check if subject folder exists;
% If it doesn't, create it and move on. p.run_num is 1.
% If it does, find the most recent run number. p.run_num is now most_recent_run_num + 1;

cd(dirs.data_dir)

if exist(p.subj_ID,'dir') == 0
    
    mkdir(p.subj_ID)
    p.run_num = 1;
    
else
    
    data_files = dir([p.subj_ID '/' dirs.save_filename_template]);
    
    if isempty(data_files)
        
        p.run_num = 1;
        
    else
        
        run_nums = zeros(1, length(data_files));
        for i = 1:length(data_files)
            run_match = regexp(data_files(i).name, 'Run(\d+)', 'tokens');
            if ~isempty(run_match)
                run_nums(i) = str2double(run_match{1}{1});
            end
        end
        p.run_num = max(run_nums) + 1;
        
    end
    
end
cd(dirs.script_dir)

disp(['Entering Run ' num2str(p.run_num)]);
if p.training
    save_filename = ['SD_Noise_Pilot_Training_S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '.mat'];
else
    save_filename = ['SD_Noise_Pilot_S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '.mat'];
end

%% Initialize experiment

init_experiment

%% Enter experiment

experiment

%% Check frame timing

ps = plotSettings();
plot_frame_timing_stats(exe_timing, p, dirs, t, ps);

%% Save run info

loading_text = 'Saving run info...';
screen_text_boundary = Screen('TextBounds', w.window, loading_text);
screen_text_patch = CenterRectOnPoint(screen_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, screen_text_patch(1),  screen_text_patch(2), w.white);
Screen('Flip', w.window);

run_info.behav_data = behav_data;
run_info.p = p;
run_info.t = t;
run_info.w = w;
run_info.frames = frames;
% run_info.pres_timing = pres_timing;
% run_info.exe_timing = exe_timing;

cd(dirs.data_dir)

if ~exist(p.subj_ID,'dir')
    mkdir(p.subj_ID)
    disp(['/' p.subj_ID ' created.'])
end

save([p.subj_ID '/' save_filename],'run_info','-mat','-v7.3');

disp(['Run #' num2str(p.run_num) ' saved!']);
disp(['Next run #: ' num2str(p.run_num +1)]);

cd(dirs.script_dir)

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

end