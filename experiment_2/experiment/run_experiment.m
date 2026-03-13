%%% run_experiment

%{

Written by Luis D. Ramirez & Natalya Phillips
UCSD
lur003@ucsd.edu

%}


function run_info = run_experiment(p,w,dirs)

t_session_start = GetSecs;

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
    save_filename = ['SD_Noise_Exp2_Training_S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '.mat'];
elseif p.calibration
    save_filename = ['SD_Noise_Exp2_Calibration_S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '.mat'];
else
    save_filename = ['SD_Noise_Exp2_S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '.mat'];
end

%%  Initialize and enter experiment

init_experiment

experiment

%% Check frame timing

plot_frame_timing_stats(exe_timing, p, dirs, t);

%% Save run info

loading_text = 'Saving run info...';
screen_text_boundary = Screen('TextBounds', w.window, loading_text);
screen_text_patch = CenterRectOnPoint(screen_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, screen_text_patch(1),  screen_text_patch(2), w.white);
Screen('Flip', w.window);

t.session_start_time = t_session_start;
t.session_end_time = GetSecs;
t.session_dur = (t.session_end_time - t.session_start_time) / 60; % min

if p.disp_on
    disp(['Actual total session duration: ' num2str(round(t.session_dur, 2)) ' minutes']);
end

run_info.behav_data = behav_data;
run_info.p = p;
run_info.t = t;
run_info.w = w;
run_info.frames = frames;
run_info.exe_timing = exe_timing;

cd(dirs.data_dir)

if ~exist(p.subj_ID,'dir')
    mkdir(p.subj_ID)
    disp(['/' p.subj_ID ' created.'])
end

save([p.subj_ID '/' save_filename],'run_info','-mat','-v7.3');

disp(['Run #' num2str(p.run_num) ' saved!']);
disp(['Next run #: ' num2str(p.run_num +1)]);

cd(dirs.script_dir)

end