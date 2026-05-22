function packed = generateClusterSample(tbl_trials, bootstrap_plan_row)
% generateClusterSample  Subject-cluster resample: one row of indices into unique subjects.
%
% tbl_trials must have subject_id. bootstrap_plan_row is 1 x N vector of integers in 1:n_subj
% indicating which subject (by stable unique order) is duplicated at each of N slots.
%
% Returns struct with numeric arrays for optimization (no table inside parfor hot path).

subj_list = unique(tbl_trials.subject_id, 'stable');
n_subj = numel(subj_list);

% Row index sets per subject (column cell)
row_sets = cell(n_subj, 1);
for s = 1:n_subj
    row_sets{s} = findSubjectRows(tbl_trials.subject_id, subj_list(s));
end

pool_idx = [];
boot_subject_id = [];
for j = 1:numel(bootstrap_plan_row)
    s = bootstrap_plan_row(j);
    pool_idx = [pool_idx; row_sets{s}]; %#ok<AGROW>
    boot_subject_id = [boot_subject_id; repmat(j, numel(row_sets{s}), 1)]; %#ok<AGROW>
end

packed.delta_theta = tbl_trials.delta_theta(pool_idx);
packed.x_probe = tbl_trials.x_probe(pool_idx);
packed.response = tbl_trials.response(pool_idx);
packed.subject_id = boot_subject_id;
packed.cond_manipulation = tbl_trials.cond_manipulation(pool_idx);
packed.cond_prev = tbl_trials.cond_prev(pool_idx);
packed.cond_curr = tbl_trials.cond_curr(pool_idx);

n = numel(packed.delta_theta);
packed.manipulation = ones(n, 1);
packed.manipulation(packed.cond_manipulation == 'precision') = 2;
end

function idx = findSubjectRows(subject_col, sid)
idx = find(subject_col == sid);
end
