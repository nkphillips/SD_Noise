function [sd_ci_subject, sd_subject] = bootstrapSerialDependenceSubjects(sd, num, bootstrap, p)
% bootstrapSerialDependenceSubjects
%
% Compute subject-bootstrap confidence intervals for serial-dependence
% parameters from already-estimated subject-level DoG parameters. The
% subject is the resampling unit; each bootstrap sample draws num.subjs
% subjects with replacement and averages the parameter over sampled subjects.

    if nargin < 4
        p = struct();
    end
    if ~isfield(bootstrap, 'B_subject_sd') || isempty(bootstrap.B_subject_sd)
        if isfield(bootstrap, 'B') && ~isempty(bootstrap.B)
            bootstrap.B_subject_sd = bootstrap.B;
        else
            bootstrap.B_subject_sd = 5;
        end
    end
    if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
        bootstrap.ci = [2.5, 97.5];
    end

    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
        num_params = 3;
    else
        num_params = size(sd.all.params_est, 4);
    end

    sz = [num.levels, num.levels, num.conds, num_params];
    sd_ci_subject.lo = nan(sz);
    sd_ci_subject.hi = nan(sz);
    sd_ci_subject.B = bootstrap.B_subject_sd;
    sd_ci_subject.ci = bootstrap.ci;
    sd_ci_subject.method = 'subject_percentile';

    sd_subject.mean = nan(sz);
    sd_subject.sem = nan(sz);
    sd_subject.n = nan(sz);
    sd_subject.method = 'subject_level_params';

    B = bootstrap.B_subject_sd;
    prc = bootstrap.ci;

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                for param = 1:num_params
                    subj_vals = nan(num.subjs, 1);
                    for subj = 1:num.subjs
                        subj_vals(subj) = sd.ind(subj).params_est(prev_lvl, curr_lvl, cond, param);
                    end

                    finite_vals = subj_vals(isfinite(subj_vals));
                    n_finite = numel(finite_vals);
                    sd_subject.n(prev_lvl, curr_lvl, cond, param) = n_finite;

                    if n_finite == 0
                        continue
                    end

                    sd_subject.mean(prev_lvl, curr_lvl, cond, param) = mean(finite_vals, 'omitnan');
                    if n_finite > 1
                        sd_subject.sem(prev_lvl, curr_lvl, cond, param) = ...
                            std(finite_vals, 'omitnan') / sqrt(n_finite);
                    else
                        sd_subject.sem(prev_lvl, curr_lvl, cond, param) = NaN;
                    end

                    boot_means = nan(B, 1);
                    for b = 1:B
                        sample_idx = randi(num.subjs, num.subjs, 1);
                        sample_vals = subj_vals(sample_idx);
                        boot_means(b) = mean(sample_vals, 'omitnan');
                    end

                    boot_means = boot_means(isfinite(boot_means));
                    if ~isempty(boot_means)
                        q = prctile(boot_means, prc);
                        sd_ci_subject.lo(prev_lvl, curr_lvl, cond, param) = q(1);
                        sd_ci_subject.hi(prev_lvl, curr_lvl, cond, param) = q(2);
                    end
                end
            end
        end
    end
end
