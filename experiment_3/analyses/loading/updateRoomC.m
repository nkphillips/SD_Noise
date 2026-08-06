% updateRoomC
% For each subject, update run_info.toggles.which_setup to a value of 2 from 1
%   - Load each subject's run_info
%       - hints: dir(), load(), save()
%   - Directly update run_info.toggles.which_setup
%   - Save run_info into the same 'path/filename'
%
% see:
% /loading/load_matfiles.m
% /experiment_1/experiment/run_experiment.m (this is where run_info is defined and saved (toward the bottom of the script))



%%% Load matfiles from experimental runs


%%% updateRoomC
% For each subject, update run_info.toggles.which_setup to 2

clear all; clc;

data_dir = '../../data/';    
cd(data_dir)

p.subj_IDs = {'001', '002' '003', '004' '006', '007', '008', '009', '010' '011', '013', '014', '015'};
which_setup = '3329C_ASUS';

for subj = 1:numel(p.subj_IDs)
    
    curr_subj = p.subj_IDs{subj};
    disp(['Fixing ' curr_subj '...'])

    %% Find subject data

    data_files = dir(curr_subj);
    data_file_names = {data_files.name};
    data_file_names(~contains(data_file_names, which_setup) | contains(data_file_names, {'Training', 'staircase'})) = [];


    %% Load + update + save each run_info


    cd(curr_subj)

    for n_file = 1:numel(data_file_names)

        load(data_file_names{n_file}, 'run_info'); % loads run_info struct

        run_info.toggles.which_setup = 2;   % Update which_setup
        
        save(data_file_names{n_file}, 'run_info', '-mat', '-v7.3'); % SAVE BACK

    end

    cd('..')

    disp([curr_subj ' is fixed!'])

end

