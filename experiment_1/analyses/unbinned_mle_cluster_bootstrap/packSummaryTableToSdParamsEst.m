function params_est = packSummaryTableToSdParamsEst(summary_table)
% packSummaryTableToSdParamsEst  Map results.summary_table to legacy sd.all.params_est shape.
%
% Output params_est is 3 x 3 x 2 x 4 with indices (prev_lvl, curr_lvl, manipulation, param):
%   (:,:,:,1) = A (pooled full-data point estimate)
%   (:,:,:,2) = w in 1/deg (pooled full-data point estimate)
%   (:,:,:,3) = baseline beta for calcDoG (pooled full-data point estimate)
%   (:,:,:,4) = psychometric sigma (pooled full-data point estimate)
%
% Manipulation dim: 1 = contrast, 2 = precision.

    params_est = nan(3, 3, 2, 4);

    for r = 1:height(summary_table)
        prev = summary_table.cond_prev(r);
        curr = summary_table.cond_curr(r);
        if strcmpi(char(summary_table.cond_manipulation(r)), 'contrast')
            m = 1;
        else
            m = 2;
        end

        params_est(prev, curr, m, 1) = summary_table.A_point(r);
        params_est(prev, curr, m, 2) = summary_table.w_point(r);
        params_est(prev, curr, m, 3) = summary_table.beta_point(r);
        params_est(prev, curr, m, 4) = summary_table.sigma_point(r);
    end

end
