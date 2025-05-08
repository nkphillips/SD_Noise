%%% Load matfiles from experimental runs

%% Pre-define arrays

array_size = cell(1, num.subjs);

all_runs = array_size; % will store all runs from each subject

num.trials_per_block = array_size;
num.runs = nan(1, num.subjs);

for subj = 1:num.subjs

    cd(data_dir)

    %% Find subject data

    curr_subj = subj_IDs{subj};

    data_files = dir(curr_subj);
    data_file_names = {data_files.name};

    data_file_names(~contains(data_file_names, which_setup) | contains(data_file_names, {'Training', 'staircase'})) = [];

    %% Load data

    num.runs(subj) = numel(data_file_names);

    cd(curr_subj)

    for n_file = 1:num.runs(subj)

        load(data_file_names{n_file}); % loads run_info struct
        disp([data_file_names{n_file} ' loaded']);

        all_runs{subj} = [all_runs{subj}, run_info];

        if n_file == 1 && subj == 1

            % Pull experiment info
            % Info that is identical no matter the subject

            num.blocks = all_runs{subj}(n_file).p.num_blocks;
            num.trials_per_block = all_runs{subj}(n_file).p.num_trials_per_block;

        end

        clear run_info

    end

    % Experimental info that varies with subjects
    num.runs(subj) = numel(all_runs{subj});

    cd('..')

end

cd(script_dir);
