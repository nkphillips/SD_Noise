function perf = printSuperSubjectPerformance(delta_theta_windows, p)
% printSuperSubjectPerformance
% Prints super-subject performance (Percent Correct and Percent CCW)
% averages and SEM for each major condition (Contrast, Precision) and for
% each prev->curr level combination within those conditions.
%
% Inputs
%   delta_theta_windows: struct with field .all containing cell arrays
%       responses{prev, curr, cond, window} (binary 0/1; 1=CW)
%       probe_offsets{prev, curr, cond, window} (signed degrees)
%   p: struct with fields cond_names, contrast, precision
%
% Output
%   perf: struct with aggregated counts and computed metrics

    responses_cell = delta_theta_windows.all.responses;
    probe_offsets_cell = delta_theta_windows.all.probe_offsets;

    num_levels = size(responses_cell, 1);
    num_conds = size(responses_cell, 3);
    num_windows = size(responses_cell, 4);

    % Storage
    perf = struct();
    perf.num_levels = num_levels;
    perf.num_conds = num_conds;
    perf.num_windows = num_windows;

    % Helper to accumulate counts
    function [n_trials, n_correct, n_ccw] = accumulate_counts(prev_lvl, curr_lvl, cond)
        n_trials = 0; n_correct = 0; n_ccw = 0;
        for iw = 1:num_windows
            cw_resp = responses_cell{prev_lvl, curr_lvl, cond, iw};
            probe_offs = probe_offsets_cell{prev_lvl, curr_lvl, cond, iw};
            if ~isempty(cw_resp) && ~isempty(probe_offs)
                cw_resp = cw_resp(:);
                probe_offs = probe_offs(:);
                n = min(numel(cw_resp), numel(probe_offs));
                if n == 0, continue; end
                cw_resp = cw_resp(1:n);
                probe_offs = probe_offs(1:n);
                is_correct = (cw_resp == (probe_offs > 0));
                n_trials = n_trials + n;
                n_correct = n_correct + sum(is_correct, 'omitnan');
                n_ccw = n_ccw + sum(1 - cw_resp, 'omitnan');
            end
        end
    end

    % Print header
    disp(' ');
    disp('=== SUPER SUBJECT PERFORMANCE (Overall and By Level) ===');

    for cond = 1:num_conds
        % Overall within this major condition (collapse over all prev/curr)
        total_trials = 0; total_correct = 0; total_ccw = 0;
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                [n_t, n_c, n_cwccw] = accumulate_counts(prev_lvl, curr_lvl, cond);
                total_trials = total_trials + n_t;
                total_correct = total_correct + n_c;
                total_ccw = total_ccw + n_cwccw;
            end
        end

        if total_trials > 0
            p_correct = total_correct / total_trials;
            se_correct = sqrt(max(p_correct * (1 - p_correct), 0) / total_trials) * 100;
            p_ccw = total_ccw / total_trials;
            se_ccw = sqrt(max(p_ccw * (1 - p_ccw), 0) / total_trials) * 100;
            pc_overall = 100 * p_correct;
            pccw_overall = 100 * p_ccw;
        else
            pc_overall = NaN; se_correct = NaN; pccw_overall = NaN; se_ccw = NaN;
        end

        cond_name = p.cond_names{cond};
        fprintf('\n%s (Overall):\n', cond_name);
        fprintf('  %% Correct: %.2f ± %.2f\n', pc_overall, se_correct);
        fprintf('  %% CCW    : %.2f ± %.2f\n', pccw_overall, se_ccw);

        % By level combination
        fprintf('%s (By previous -> current level):\n', cond_name);
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                [n_t, n_c, n_cwccw] = accumulate_counts(prev_lvl, curr_lvl, cond);
                if n_t > 0
                    p_corr = n_c / n_t; se_corr = sqrt(max(p_corr * (1 - p_corr), 0) / n_t) * 100; pc = 100 * p_corr;
                    p_ccw = n_cwccw / n_t; se_ccw = sqrt(max(p_ccw * (1 - p_ccw), 0) / n_t) * 100; pccw = 100 * p_ccw;
                else
                    pc = NaN; se_corr = NaN; pccw = NaN; se_ccw = NaN;
                end
                if cond == 1
                    prev_label = p.contrast{prev_lvl}; curr_label = p.contrast{curr_lvl};
                else
                    prev_label = p.precision{prev_lvl}; curr_label = p.precision{curr_lvl};
                end
                fprintf('  %s -> %s | %%Correct: %.2f ± %.2f, %%CCW: %.2f ± %.2f (N=%d)\n', ...
                    prev_label, curr_label, pc, se_corr, pccw, se_ccw, n_t);
                % Save into struct
                perf.by_level(prev_lvl, curr_lvl, cond).N = n_t;
                perf.by_level(prev_lvl, curr_lvl, cond).pc = pc;
                perf.by_level(prev_lvl, curr_lvl, cond).pc_sem = se_corr;
                perf.by_level(prev_lvl, curr_lvl, cond).pccw = pccw;
                perf.by_level(prev_lvl, curr_lvl, cond).pccw_sem = se_ccw;
            end
        end

        perf.overall(cond).N = total_trials;
        perf.overall(cond).pc = pc_overall;
        perf.overall(cond).pc_sem = se_correct;
        perf.overall(cond).pccw = pccw_overall;
        perf.overall(cond).pccw_sem = se_ccw;
    end

    disp('===============================================');
end


