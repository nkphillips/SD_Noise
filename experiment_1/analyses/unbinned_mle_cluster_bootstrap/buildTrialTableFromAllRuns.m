function tbl_trials = buildTrialTableFromAllRuns(all_runs, p, num, n_back)
% buildTrialTableFromAllRuns  Trial-level table for unbinned MLE (no delta-theta windows).

    subj_id_col = [];
    cond_manip_col = {};
    cond_prev_col = [];
    cond_curr_col = [];
    delta_col = [];
    x_probe_col = [];
    resp_col = [];

    for subj = 1:num.subjs

        for n_run = 1:num.runs(subj)

            subj_p = all_runs{subj}(n_run).p;
            behav_data = all_runs{subj}(n_run).behav_data;

            delta_thetas = calcOrientationDiff(squeeze(subj_p.trial_events(1:end-n_back, 1, :)), squeeze(subj_p.trial_events(n_back+1:end, 1, :)));
            probe_offsets = calcOrientationDiff(subj_p.trial_events(:, 2, :), subj_p.trial_events(:, 1, :));

            for cond = 1:num.conds

                curr_cond_blocks = find(subj_p.cond_order == cond);
                if cond == 1
                    manip_label = 'contrast';
                else
                    manip_label = 'precision';
                end

                for n_block = 1:length(curr_cond_blocks)

                    curr_delta_thetas = delta_thetas(:, curr_cond_blocks(n_block));
                    curr_lvls = subj_p.trial_events(:, 3, curr_cond_blocks(n_block));
                    curr_probe_offsets = probe_offsets(n_back+1:end, curr_cond_blocks(n_block));
                    curr_CW_response = behav_data.response(n_back+1:end, curr_cond_blocks(n_block)) == 2;

                    for prev_lvl = 1:num.levels
                        for curr_lvl = 1:num.levels

                            curr_lvl_pair_indx = curr_lvls(1:end-n_back) == prev_lvl & curr_lvls(n_back+1:end) == curr_lvl;

                            if ~any(curr_lvl_pair_indx)
                                continue
                            end

                            n_t = sum(curr_lvl_pair_indx);
                            subj_id_col = [subj_id_col; repmat(subj, n_t, 1)]; %#ok<AGROW>
                            cond_manip_col = [cond_manip_col; repmat({manip_label}, n_t, 1)]; %#ok<AGROW>
                            cond_prev_col = [cond_prev_col; repmat(prev_lvl, n_t, 1)]; %#ok<AGROW>
                            cond_curr_col = [cond_curr_col; repmat(curr_lvl, n_t, 1)]; %#ok<AGROW>
                            delta_col = [delta_col; curr_delta_thetas(curr_lvl_pair_indx)]; %#ok<AGROW>
                            x_probe_col = [x_probe_col; curr_probe_offsets(curr_lvl_pair_indx)]; %#ok<AGROW>
                            resp_col = [resp_col; double(curr_CW_response(curr_lvl_pair_indx))]; %#ok<AGROW>

                        end
                    end

                end

            end

        end

    end

    tbl_trials = table(subj_id_col, categorical(cond_manip_col), cond_prev_col, cond_curr_col, ...
        delta_col, x_probe_col, resp_col, ...
        'VariableNames', {'subject_id', 'cond_manipulation', 'cond_prev', 'cond_curr', ...
        'delta_theta', 'x_probe', 'response'});

end
