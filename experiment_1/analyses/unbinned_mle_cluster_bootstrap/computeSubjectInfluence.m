function influence = computeSubjectInfluence(results)
% computeSubjectInfluence  Per-subject leverage on each per-cell parameter,
% computed from the existing leave-one-subject-out jackknife.
%
% influence_i,c,k = theta_full(c,k) - theta_jk_i(c,k)
%
% Where theta_full is the pooled-data point estimate (results.overlay.params_point)
% and theta_jk_i is the leave-one-out estimate excluding subject i
% (results.jackknife.params).
%
% Returns:
%   influence.values             [n_subj x num_conds x 4]  signed influence
%   influence.values_abs         [n_subj x num_conds x 4]  |influence|
%   influence.median_abs         [n_subj x 4]              median |influence| across 18 cells
%   influence.mean_abs           [n_subj x 4]              mean   |influence|
%   influence.max_abs            [n_subj x 4]              max    |influence|
%   influence.l2                 [n_subj x 4]              sqrt(sum(influence^2))
%   influence.flagged            [n_subj x 4]  logical     median_abs above (median + 2*IQR) threshold
%   influence.flag_threshold     [1 x 4]                   threshold value used for flagging
%   influence.param_labels       {'A','w','sigma','beta'}  for column reference

    if ~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params')
        error('computeSubjectInfluence:noJackknife', ...
            'Requires results.jackknife.params; rerun with compute_jackknife = true.');
    end

    pp = results.overlay.params_point;       % num_conds x 4
    jk = results.jackknife.params;            % n_subj x num_conds x 4
    [n_subj, num_conds, n_params] = size(jk);

    values = nan(n_subj, num_conds, n_params);
    for k = 1:n_params
        for i = 1:n_subj
            values(i, :, k) = pp(:, k)' - squeeze(jk(i, :, k));
        end
    end

    influence = struct();
    influence.values = values;
    influence.values_abs = abs(values);
    influence.median_abs = squeeze(median(abs(values), 2, 'omitnan'));
    influence.mean_abs   = squeeze(mean(abs(values), 2, 'omitnan'));
    influence.max_abs    = squeeze(max(abs(values), [], 2));
    influence.l2         = squeeze(sqrt(sum(values.^2, 2, 'omitnan')));

    flagged = false(n_subj, n_params);
    flag_threshold = nan(1, n_params);
    for k = 1:n_params
        m = median(influence.median_abs(:, k), 'omitnan');
        iqr_val = iqr(influence.median_abs(:, k));
        flag_threshold(k) = m + 2 * iqr_val;
        flagged(:, k) = influence.median_abs(:, k) > flag_threshold(k);
    end
    influence.flagged = flagged;
    influence.flag_threshold = flag_threshold;
    influence.param_labels = {'A', 'w', 'sigma', 'beta'};
end
