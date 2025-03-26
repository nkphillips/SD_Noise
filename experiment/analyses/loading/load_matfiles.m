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

    data_file_names(~contains(data_file_names, which_setup) | contains(data_file_names, 'Training')) = [];

    %% Load data

    num.runs(subj) = numel(data_file_names);

    cd(curr_subj)

    for n_file = 1:num.runs(subj)

        load(data_file_names{n_file}); % loads run_info struct

        all_runs{subj} = [all_runs{subj}, run_info];

        num.trials_per_block{subj}(n_file,:) = run_info.behav_data.num_trials_per_block;

        if n_file == 1 && subj == 1

            % Pull experiment info
            % Info that is identical no matter the subject

            num.test_sfs = all_runs{subj}(n_file).p.num_test_sfs;
            num.blocks = all_runs{subj}(n_file).p.num_blocks;

            num.min_trials_per_block = all_runs{subj}(n_file).p.num_trials_per_block;

            num.min_trials_per_cond_per_run = num.min_trials_per_block * (num.blocks/num.conds);
            num.min_trials_per_test_sf_per_cond_per_run = num.min_trials_per_block/num.test_sfs;

            num.min_total_trials_per_cond = num.min_trials_per_cond_per_run * num.runs(subj);
            num.min_total_trials = num.min_total_trials_per_cond * num.conds;

            num.trial_per_test_sf_per_cond = num.min_trials_per_test_sf_per_cond_per_run * num.runs(subj);

            test_sfs = all_runs{subj}(n_file).p.test_sfs;
            reference_sf = all_runs{subj}(n_file).p.reference_sf;
            adaptor_sfs = all_runs{subj}(n_file).p.adaptor_sfs;
            sf_ratios = all_runs{subj}(n_file).p.sf_ratios;

            delta_sfs = nan(num.test_sfs);
            for i = 1:num.test_sfs
                for j = 1:num.test_sfs

                    delta_sfs(i,j) = log2(test_sfs(i)/test_sfs(j));

                end
            end
        end

        clear run_info

    end

    % Experimental info that varies with subjects
    num.runs(subj) = numel(all_runs{subj});
    num.trials_per_cond_per_run{subj} = num.trials_per_block{subj} * (num.blocks/num.conds);
    num.total_trials_per_cond{subj} = num.trials_per_cond_per_run{subj} * num.runs(subj);
    num.total_trials{subj} = num.total_trials_per_cond{subj} * num.conds;

    cd('..')

end

cd(script_dir);
