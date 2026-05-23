function perf = computeDeltaThetaPerformanceFromTrialTable(tbl_trials, delta_theta_centers, delta_theta_width)
% computeDeltaThetaPerformanceFromTrialTable  Legacy plotPerformance inputs from trial tables.
%
% Computes percent correct and percent CCW in sliding delta-theta windows for
% the 3 x 3 previous/current level cells of each manipulation. The output fields
% match the arrays consumed by plotting/plotPerformance.m.

    if nargin < 3 || isempty(delta_theta_width)
        delta_theta_width = 32;
    end

    num_levels = 3;
    num_conds = 2;
    num_windows = numel(delta_theta_centers);

    perf.performance = nan(num_levels, num_levels, num_conds, num_windows);
    perf.pCW = nan(num_levels, num_levels, num_conds, num_windows);
    perf.n_trials = zeros(num_levels, num_levels, num_conds, num_windows);
    perf.delta_theta_centers = delta_theta_centers(:)';
    perf.delta_theta_width = delta_theta_width;

    for prev_lvl = 1:num_levels
        for curr_lvl = 1:num_levels
            for cond = 1:num_conds
                if cond == 1
                    manip = 'contrast';
                else
                    manip = 'precision';
                end

                base_mask = tbl_trials.cond_manipulation == manip & ...
                    tbl_trials.cond_prev == prev_lvl & ...
                    tbl_trials.cond_curr == curr_lvl & ...
                    (tbl_trials.response == 0 | tbl_trials.response == 1) & ...
                    isfinite(tbl_trials.delta_theta) & ...
                    isfinite(tbl_trials.x_probe);

                for i_window = 1:num_windows
                    left_edge = delta_theta_centers(i_window) - delta_theta_width / 2;
                    right_edge = delta_theta_centers(i_window) + delta_theta_width / 2;

                    if left_edge < -90
                        in_window = tbl_trials.delta_theta >= left_edge + 180 | ...
                            tbl_trials.delta_theta <= right_edge;
                    elseif right_edge > 90
                        in_window = tbl_trials.delta_theta <= right_edge - 180 | ...
                            tbl_trials.delta_theta >= left_edge;
                    else
                        in_window = tbl_trials.delta_theta >= left_edge & ...
                            tbl_trials.delta_theta <= right_edge;
                    end

                    mask = base_mask & in_window;
                    n = sum(mask);
                    perf.n_trials(prev_lvl, curr_lvl, cond, i_window) = n;

                    if n > 0
                        resp = tbl_trials.response(mask);
                        x_probe = tbl_trials.x_probe(mask);
                        is_correct = resp == (x_probe > 0);
                        perf.performance(prev_lvl, curr_lvl, cond, i_window) = ...
                            100 * mean(is_correct, 'omitnan');
                        perf.pCW(prev_lvl, curr_lvl, cond, i_window) = ...
                            100 * mean(1 - resp, 'omitnan');
                    end
                end
            end
        end
    end
end
