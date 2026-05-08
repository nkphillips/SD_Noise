load('experiment_2/data/000/SD_Noise_Exp2_Calibration_S000_Run1_3329C_ASUS.mat', 'run_info');
r000 = run_info;
load('experiment_2/data/777/SD_Noise_Exp2_Calibration_S777_Run1_3329C_ASUS.mat', 'run_info');
r777 = run_info;

disp('Are responses the same?');
isequal(r000.behav_data.response, r777.behav_data.response)

disp('Are corrects the same?');
isequal(r000.behav_data.correct, r777.behav_data.correct)

disp('Are trial events the same?');
isequal(r000.p.trial_events, r777.p.trial_events)

