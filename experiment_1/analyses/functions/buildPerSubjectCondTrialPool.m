function trial_pool = buildPerSubjectCondTrialPool(all_runs, p, num, n_back) %#ok<INUSL>
% buildPerSubjectCondTrialPool  Per-(prev,curr,cond) trial vectors WITHOUT Δθ windowing.
%
% Returns a struct with .ind(subj) and .all per-condition cell arrays of trial-level
% data, each trial counted exactly ONCE. This is the data used for trial-level NLL
% fits (S&S-style joint MLE of A, w, sigma over all trials in a (manip, prev, curr)
% cell).
%
% Compare to makeDeltaThetaWindows: that function bins each trial into ~32 overlapping
% Δθ windows; vertcatting across windows duplicates each trial ~32×, which inflates
% likelihood-based fits. Use this helper for any NLL or BCa step that pools trials
% directly; use makeDeltaThetaWindows.m for response-bias/SSE paths that genuinely
% need per-window summaries.
%
% Output:
%   trial_pool.ind(subj).delta_thetas   {prev_lvl, curr_lvl, cond}  Nx1 vector
%   trial_pool.ind(subj).probe_offsets  {prev_lvl, curr_lvl, cond}  Nx1 vector
%   trial_pool.ind(subj).responses      {prev_lvl, curr_lvl, cond}  Nx1 vector (0 = CCW, 1 = CW)
%   trial_pool.all.delta_thetas         {prev_lvl, curr_lvl, cond}  super-subject pool

    fields = {'delta_thetas', 'probe_offsets', 'responses'};

    template_ind = struct();
    for i_field = 1:numel(fields)
        template_ind.(fields{i_field}) = cell(num.levels, num.levels, num.conds);
    end

    trial_pool.ind = repmat(template_ind, num.subjs, 1);
    trial_pool.all = template_ind;

    for subj = 1:num.subjs
        for n_run = 1:num.runs(subj)

            subj_p = all_runs{subj}(n_run).p;
            behav_data = all_runs{subj}(n_run).behav_data;

            delta_thetas = calcOrientationDiff( ...
                squeeze(subj_p.trial_events(1:end-n_back, 1, :)), ...
                squeeze(subj_p.trial_events(n_back+1:end, 1, :)));
            probe_offsets = calcOrientationDiff( ...
                subj_p.trial_events(:, 2, :), subj_p.trial_events(:, 1, :));

            for cond = 1:num.conds

                curr_cond_blocks = find(subj_p.cond_order == cond);

                for n_block = 1:length(curr_cond_blocks)

                    blk = curr_cond_blocks(n_block);
                    curr_delta_thetas  = delta_thetas(:, blk);
                    curr_lvls          = subj_p.trial_events(:, 3, blk);
                    curr_probe_offsets = probe_offsets(n_back+1:end, blk);
                    curr_CW_response   = behav_data.response(n_back+1:end, blk) == 2;

                    for prev_lvl = 1:num.levels
                        for curr_lvl = 1:num.levels

                            sel = curr_lvls(1:end-n_back) == prev_lvl & ...
                                  curr_lvls(n_back+1:end) == curr_lvl;

                            if ~any(sel)
                                continue
                            end

                            dt = curr_delta_thetas(sel);
                            po = curr_probe_offsets(sel);
                            rv = double(curr_CW_response(sel));

                            trial_pool.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond}; dt(:)];
                            trial_pool.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond}; po(:)];
                            trial_pool.ind(subj).responses{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.ind(subj).responses{prev_lvl, curr_lvl, cond}; rv(:)];

                            trial_pool.all.delta_thetas{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.all.delta_thetas{prev_lvl, curr_lvl, cond}; dt(:)];
                            trial_pool.all.probe_offsets{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.all.probe_offsets{prev_lvl, curr_lvl, cond}; po(:)];
                            trial_pool.all.responses{prev_lvl, curr_lvl, cond} = ...
                                [trial_pool.all.responses{prev_lvl, curr_lvl, cond}; rv(:)];

                        end
                    end

                end

            end

        end
    end

end
