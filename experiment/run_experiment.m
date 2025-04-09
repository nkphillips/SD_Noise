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
        
        data_file_names = sort({data_files.name});
        most_recent_file = data_file_names{end};
        
        file_strs = strsplit(most_recent_file,'_');
        if ~p.training
            most_recent_run_num = str2double(file_strs{5}(end));
        else
            most_recent_run_num = str2double(file_strs{6}(end));
        end
        
        p.run_num = most_recent_run_num + 1;
        
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

%% Check event timing

%{

% Set figure path and check existence
figure_path = dirs.logs_dir;
if ~exist(figure_path,'dir')
    mkdir(figure_path)
    disp([figure_path ' created.'])
end

% Open figure
fg  = figure('Color', [1 1 1]);
set(0, 'CurrentFigure', fg)
figure_name = ['S' p.subj_ID ' Frame timing Run ' num2str(p.run_num) ' ' p.which_setup];

% Plot frame timing
subplot(1,2,1)
x = cat(1,exe_timing.test_Missed{:});
plot(x'); hold on; line([min(xlim) max(xlim)], [0 0],'Color', 'k', 'LineStyle', '-')
test_frames_missed = sum(x(:) > 0);
title(['Test frames missed: ' num2str(test_frames_missed) ' (' num2str(round(mean(x(:) > 0) * 100)) '%)'])
box off; axis square;
xlabel('Frame #')
ylabel('Deadline offset (s)')

subplot(1,2,2)
x = cat(1,exe_timing.top_up_Missed{:});
plot(x'); hold on; line([min(xlim) max(xlim)], [0 0], 'Color', 'k', 'LineStyle', '-')
top_up_frames_missed = sum(x(:) > 0);
title(['Top up frames missed: ' num2str(top_up_frames_missed) ' (' num2str(round(mean(x(:) > 0) * 100)) '%)'])
box off; axis square;
xlabel('Frame #')
ylabel('Deadline offset (s)')

saveas(gcf,[figure_path '/' figure_name '.pdf']);
close(fg)

%}

%% Save run info

loading_text = 'Saving run info...';
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
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
    
    loading_text = 'Saving textures...';
    loading_text_boundary = Screen('TextBounds', w.window, loading_text);
    loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);
    
    Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
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
    
    % Check if subject data folder exists
    if ~exist(p.subj_ID, 'dir')
        mkdir(p.subj_ID)
    end
    
    cd(p.subj_ID)
    if ~exist(textures_filename, "file")
        
        % If textures don't exist, save them
        disp('Saving textures matfile...')
        save(textures_filename,'textures','-mat','-v7.3');
        
    else
        
        % If they do exist, manually decide whether to overwrite
        if ~p.simulate_response
            overwrite_file = [];
        else
            overwrite_file = 0;
        end
                    
        loading_text = [textures_filename ' exists. Do you want to overwrite it? (1 - yes, 0 - no):'];
        loading_text_boundary = Screen('TextBounds', w.window, loading_text);
        loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);
        
        while isempty(overwrite_file)

            Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
            Screen('Flip', w.window);
            
            [key_pressed, first_press] = KbQueueCheck(p.device_number);
            which_press = find(first_press);
            
            if key_pressed && which_press(1) == KbName('1!')
                overwrite_file = 1;
            elseif key_pressed && which_press(1) == KbName('0)')
                overwrite_file = 0;
            end
            
        end
        
        if overwrite_file
            save(textures_filename,'textures','-mat','-v7.3');
        end
        
    end
    
    loading_text = 'Textures saved!';
    loading_text_boundary = Screen('TextBounds', w.window, loading_text);
    loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);
    
    Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
    Screen('Flip', w.window);
    
    cd(dirs.script_dir)
    
end

end