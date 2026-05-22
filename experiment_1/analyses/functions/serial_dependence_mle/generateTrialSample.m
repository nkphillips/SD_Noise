function packed = generateTrialSample(tbl_trials, bootstrap_plan_row)
% generateTrialSample  Trial-level resample from an already-built history table.
%
% The trial table already encodes within-subject/run/block history. This helper
% resamples those rows as observations without recomputing temporal history.

    pool_idx = bootstrap_plan_row(:);

    packed.delta_theta = tbl_trials.delta_theta(pool_idx);
    packed.x_probe = tbl_trials.x_probe(pool_idx);
    packed.response = tbl_trials.response(pool_idx);
    packed.subject_id = tbl_trials.subject_id(pool_idx);
    packed.cond_manipulation = tbl_trials.cond_manipulation(pool_idx);
    packed.cond_prev = tbl_trials.cond_prev(pool_idx);
    packed.cond_curr = tbl_trials.cond_curr(pool_idx);

    n = numel(packed.delta_theta);
    packed.manipulation = ones(n, 1);
    packed.manipulation(packed.cond_manipulation == 'precision') = 2;
end
