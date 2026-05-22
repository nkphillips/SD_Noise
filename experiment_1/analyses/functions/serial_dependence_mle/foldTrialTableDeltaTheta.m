function [tbl_folded, fold_info] = foldTrialTableDeltaTheta(tbl_trials)
% foldTrialTableDeltaTheta  Fold negative delta-theta trials for odd-symmetric fits.
%
% Negative-delta trials are reflected into positive delta-theta space. The
% probe offset and binary CW response are flipped so the likelihood remains
% equivalent to fitting an odd-symmetric mu(delta_theta) curve.

    required_vars = {'delta_theta', 'x_probe', 'response'};
    for i_var = 1:numel(required_vars)
        if ~ismember(required_vars{i_var}, tbl_trials.Properties.VariableNames)
            error('foldTrialTableDeltaTheta:missingVariable', ...
                'tbl_trials is missing required variable: %s', required_vars{i_var});
        end
    end

    tbl_folded = tbl_trials;

    is_neg = tbl_folded.delta_theta < 0;
    is_valid_response = tbl_folded.response == 0 | tbl_folded.response == 1;
    flip_response = is_neg & is_valid_response;

    tbl_folded.delta_theta = abs(tbl_folded.delta_theta);
    tbl_folded.x_probe(is_neg) = -tbl_folded.x_probe(is_neg);
    tbl_folded.response(flip_response) = 1 - tbl_folded.response(flip_response);

    fold_info = struct();
    fold_info.n_trials = height(tbl_trials);
    fold_info.n_flipped = sum(is_neg);
    fold_info.n_nonnegative = sum(~is_neg);
    fold_info.n_invalid_response_on_negative_delta = sum(is_neg & ~is_valid_response);
end
