%%% init_behav_perf

behav_perf = struct('ind',[],'grp',[]);

tmp_array = cell(1, num.subjs);
behav_perf.ind = struct('mean_within_runs', tmp_array, 'sem_within_runs', tmp_array, ...
    'mean_within_test_sf', tmp_array, 'sem_within_test_sf', tmp_array, ...
    'mean_across_runs', tmp_array, 'sem_across_runs', tmp_array, ...
    'mean_across_conds', tmp_array, 'sem_across_conds', tmp_array, ...
    'mean_within_test_sf_across_runs', tmp_array, 'sem_within_test_sf_across_runs', tmp_array);

tmp_array = cell(1, 1);
behav_perf.grp = struct('mean_within_runs', tmp_array, 'sem_within_runs', tmp_array, ...
    'mean_within_test_sf', tmp_array, 'sem_within_test_sf', tmp_array, ...
    'mean_across_runs', tmp_array, 'sem_across_runs', tmp_array, ...
    'mean_across_conds', tmp_array, 'sem_across_conds', tmp_array, ...
    'mean_within_test_sf_across_runs', tmp_array, 'sem_within_test_sf_across_runs', tmp_array);

struct_fieldnames = fieldnames(behav_perf.ind);

%% Subject performance

for subj = 1:num.subjs

    % Initialize temporary matrices to store statistics
    within_run_mean = nan(num.runs(subj), num.conds);
    within_run_sem = within_run_mean;

    within_test_sf_mean = nan(num.test_sfs, num.conds, num.runs(subj));
    within_test_sf_sem = within_test_sf_mean;

    mean_within_delta_sf = nan(num.test_sfs,num.conds,num.test_sfs,num.runs(subj));
    sem_within_delta_sf = mean_within_delta_sf;

    % Get performance and sem for each run
    for n_run = 1:num.runs(subj)

        % max_num_trials = max(num.trials_per_block{subj}(n_run,:));

        for cond = 1:num.conds

            % Get block indices from current adaptation condition (ignoring hemifield)
            curr_blocks = find(all_runs{subj}(n_run).p.adapt_cond_order(1,:) == cond);

            % Pre-allocate matrix for responses and test_sf_order (b/c subjs can redo trials, the num cols varies)
            responses = nan(length(curr_blocks), num.min_trials_per_block);
            test_sf_order = nan(length(curr_blocks), num.min_trials_per_block);
            deltas = nan(length(curr_blocks), num.min_trials_per_block-1);

            % Loop through every block and pull out responses, test SF order, and delta between previous and current test SF
            for n_block = 1:length(curr_blocks)

                missed_trials = all_runs{subj}(n_run).behav_data.missed_trial_indices{curr_blocks(n_block)};

                % Get responses, test SF order, and delta SF
                curr_responses = all_runs{subj}(n_run).behav_data.response{curr_blocks(n_block)};
                curr_test_sf_order = all_runs{subj}(n_run).p.test_sf_order(curr_blocks(n_block),:);
                curr_delta = log2(curr_test_sf_order(1:end-1)./curr_test_sf_order(2:end));

                responses(n_block, 1:length(curr_responses)) = curr_responses;
                test_sf_order(n_block, 1:length(curr_test_sf_order)) = curr_test_sf_order;
                deltas(n_block, 1:length(curr_delta)) = curr_delta;

            end

            % Store performance (avg) and sem
            within_run_mean(n_run,cond) = mean(responses(:),'omitnan'); % collapses across test sfs
            within_run_sem(n_run,cond) = std(responses(:),'omitnan')/sqrt(sum(~isnan(responses(:))));

            for n_sf = 1:num.test_sfs

                test_sf_indx = test_sf_order == test_sfs(n_sf);

                curr_responses = responses(test_sf_indx);

                within_test_sf_mean(n_sf,cond,n_run) = mean(curr_responses, 'omitnan');
                within_test_sf_sem(n_sf,cond,n_run) = std(curr_responses, 'omitnan')/sqrt(length(curr_responses));

                % Performance as a function of delta SF between previous and present test SF
                curr_responses = responses(:,2:end); % only use trials after the first 
                curr_deltas = deltas(test_sf_indx(:,2:end));

                for m_sf = 1:num.test_sfs

                    delta_indx = curr_deltas == delta_sfs(m_sf, n_sf);
                    behav_perf.ind(subj).trials_per_delta_sf(m_sf,n_sf,cond) = sum(delta_indx);

                    mean_within_delta_sf(m_sf,cond,n_sf,n_run) = mean(curr_responses(delta_indx),'omitnan');
                    sem_within_delta_sf(m_sf,cond,n_sf,n_run) = std(curr_responses(delta_indx),'omitnan')/sqrt(length(curr_responses(delta_indx)));

                end

            end

        end
    end

    % Store performance within runs
    behav_perf.ind(subj).mean_within_runs = within_run_mean;
    behav_perf.ind(subj).sem_within_runs = within_run_sem;

    % Store performance within test SFs
    behav_perf.ind(subj).mean_within_test_sf = within_test_sf_mean;
    behav_perf.ind(subj).sem_within_test_sf = within_test_sf_sem;

    % Average performance across runs
    behav_perf.ind(subj).mean_across_runs = mean(within_run_mean,1);
    behav_perf.ind(subj).sem_across_runs = std(within_run_mean,[],1)./sqrt(num.runs(subj));

    % Average performance across conditions, across runs
    behav_perf.ind(subj).mean_across_conds = mean(mean(within_run_mean,2)); % Average performance (across conditions)
    behav_perf.ind(subj).sem_across_conds = std(mean(within_run_mean,2))./sqrt(num.runs(subj)); % sem (across conditions)

    % Average performance within test SFs, across runs
    behav_perf.ind(subj).mean_within_test_sf_across_runs = mean(within_test_sf_mean,3);
    behav_perf.ind(subj).sem_within_test_sf_across_runs = std(within_test_sf_mean,[],3)/sqrt(num.runs(subj));

    % Average performance within delta test SFs for ea. test SF and condition
    behav_perf.ind(subj).mean_within_delta_sf = mean(mean_within_delta_sf,4,'omitnan');
    behav_perf.ind(subj).sem_within_delta_sf = std(mean_within_delta_sf,[],4,'omitnan')./sqrt(sum(~isnan(mean_within_delta_sf),4));

end

%% Group performance

% Average performance within runs
behav_perf.grp.mean_within_runs = mean(cat(1,behav_perf.ind(:).mean_within_runs),1); % this is akin to a "super subject"
behav_perf.grp.sem_within_runs = std(cat(1,behav_perf.ind(:).mean_within_runs),[],1)./sqrt(size(cat(1,behav_perf.ind(:).mean_within_runs),1));

% Average performance within test SFs
behav_perf.grp.mean_within_test_sf = mean(cat(3,behav_perf.ind(:).mean_within_test_sf),3); % this is akin to a "super subject"
behav_perf.grp.sem_within_test_sf = std(cat(3,behav_perf.ind(:).mean_within_test_sf),[],3)/sqrt(size(cat(3,behav_perf.ind(:).mean_within_test_sf),3)); % this is akin to a "super subject"

% Average performance across subjects for each condition
subj_perf_across_runs = cat(1,behav_perf.ind(:).mean_across_runs);
behav_perf.grp.mean_across_runs = mean(subj_perf_across_runs,1);
behav_perf.grp.sem_across_runs = std(subj_perf_across_runs,[],1)/sqrt(num.subjs);

% Average performance across subjs, across conditions, across runs
subj_perf_across_conds = cat(1,behav_perf.ind(:).mean_across_conds);
behav_perf.grp.mean_across_conds = mean(subj_perf_across_conds);
behav_perf.grp.sem_across_conds = std(subj_perf_across_conds)/sqrt(num.subjs);

% Average performance across subjects for each test SF and condition
subj_perf_within_test_sf_across_runs = cat(3,behav_perf.ind(:).mean_within_test_sf_across_runs);
behav_perf.grp.mean_within_test_sf_across_runs = mean(subj_perf_within_test_sf_across_runs,3);
behav_perf.grp.sem_within_test_sf_across_runs = std(subj_perf_within_test_sf_across_runs,[],3)/sqrt(num.subjs);

% Average performance within delta test SFs for ea. test SF and condition
subj_perf_within_delta_sf = cat(4,behav_perf.ind(:).mean_within_delta_sf);
behav_perf.grp.mean_within_delta_sf = mean(subj_perf_within_delta_sf,4,'omitnan');
behav_perf.grp.sem_within_delta_sf = std(subj_perf_within_delta_sf,[],4,'omitnan')./sqrt(sum(~isnan(subj_perf_within_delta_sf),4));