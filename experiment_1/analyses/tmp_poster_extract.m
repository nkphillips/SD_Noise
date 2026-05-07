%%% tmp_poster_extract.m
% Reads unbinned_mle_cluster_bootstrap/poster_run.mat and prints formatted markdown
% suitable for pasting into POSTER_SCRIPT.md as the "Trial counts" and "Statistical
% results" subsections.

addpath('functions');
addpath('unbinned_mle_cluster_bootstrap');

S = load('unbinned_mle_cluster_bootstrap/poster_run.mat');

trial_counts = S.trial_counts;
total_trials_per_cond = S.total_trials_per_cond;
results = S.results;
num = S.num;
p = S.p;
n_back = S.n_back;

%% Trial counts
n_subj = num.subjs;
total_runs = sum(arrayfun(@(s) length(S.results.boot_plans(1,:)), 1));   % approximate

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('SUBSECTION: TRIAL COUNTS (paste into POSTER_SCRIPT.md)\n');
fprintf('============================================================\n');
fprintf('\n## Trial counts (n_back = %d, N = %d subjects)\n\n', n_back, n_subj);

per_subj_total = sum(total_trials_per_cond, 1);
fprintf('- **Total trials pooled**: %d (sum of usable previous-current pairs across the cohort).\n', sum(per_subj_total));
fprintf('- **Trials per subject**: median %.0f, range %.0f-%.0f.\n', median(per_subj_total), min(per_subj_total), max(per_subj_total));

per_cell_pooled = squeeze(sum(trial_counts, 4));   % levels x levels x conds
fprintf('- **Trials per (manipulation, prev, curr) cell, pooled across subjects**: median %.0f, range %.0f-%.0f (across the 18 cells).\n', ...
    median(per_cell_pooled(:)), min(per_cell_pooled(:)), max(per_cell_pooled(:)));

per_subj_per_cell = trial_counts(:);
per_subj_per_cell = per_subj_per_cell(per_subj_per_cell > 0);
fprintf('- **Trials per cell per subject**: median %.0f (IQR %.0f-%.0f), range %.0f-%.0f.\n', ...
    median(per_subj_per_cell), prctile(per_subj_per_cell, 25), prctile(per_subj_per_cell, 75), ...
    min(per_subj_per_cell), max(per_subj_per_cell));

%% Statistical results
fprintf('\n\n');
fprintf('============================================================\n');
fprintf('SUBSECTION: STATISTICAL RESULTS (paste into POSTER_SCRIPT.md)\n');
fprintf('============================================================\n');
fprintf('\n## Statistical results (n_back = %d, B = %d cluster-bootstrap replicates, BCa 95%% CI)\n\n', n_back, results.B);

T = results.summary_table;

fprintf('Bootstrap admission rate: median %.1f%% (range %.1f-%.1f%%) of %d replicates per cell.\n', ...
    100 * median(results.discard_summary.admitted_per_cond) / results.B, ...
    100 * min(results.discard_summary.admitted_per_cond) / results.B, ...
    100 * max(results.discard_summary.admitted_per_cond) / results.B, results.B);
fprintf('Discard reasons: failed_fit = %d, at_bound = %d (across all 18 cells x %d replicates).\n\n', ...
    results.discard_summary.by_reason.failed_fit, ...
    results.discard_summary.by_reason.at_bound, results.B);

% Cells where A CI excludes zero
A_lo = T.alpha_prctile_lo;
A_hi = T.alpha_prctile_hi;
A_med = T.alpha_median;
A_excludes_zero = (A_lo > 0) | (A_hi < 0);
A_attractive = (A_lo > 0);   % positive (attractive) and CI excludes zero
A_repulsive = (A_hi < 0);

fprintf('### Detection of serial dependence per cell\n\n');
fprintf('- **Cells with reliable attractive serial dependence** (BCa 95%% CI on A excludes zero, A > 0): **%d / 18**.\n', sum(A_attractive));
fprintf('- **Cells with reliable repulsive serial dependence** (BCa 95%% CI on A excludes zero, A < 0): **%d / 18**.\n', sum(A_repulsive));
fprintf('- **Cells where A CI includes zero (no reliable effect)**: %d / 18.\n', sum(~A_excludes_zero));

%% Tables: Contrast and Precision
fprintf('\n### Headline numbers per cell\n\n');

manip_labels = {'contrast', 'precision'};
panel_labels = {p.contrast, p.precision};

for m = 1:2
    name = manip_labels{m};
    if m == 1, mtitle = 'Contrast'; else, mtitle = 'Precision'; end
    labels = panel_labels{m};

    fprintf('\n**%s manipulation (n_back = %d)**\n\n', mtitle, n_back);
    fprintf('| Prev → Curr | A (deg) [95%% CI] | FWHM (deg) [95%% CI] | beta (deg) [95%% CI] | sigma (deg) | A excludes 0? |\n');
    fprintf('|---|---|---|---|---|---|\n');

    for prev = 1:3
        for curr = 1:3
            row_mask = T.cond_manipulation == name & T.cond_prev == prev & T.cond_curr == curr;
            r = find(row_mask, 1);
            if isempty(r), continue; end

            A   = T.alpha_median(r);
            Alo = T.alpha_prctile_lo(r);
            Ahi = T.alpha_prctile_hi(r);
            F   = unbinnedWtoFwhm(T.w_median(r));
            Flo = T.fwhm_prctile_lo(r);
            Fhi = T.fwhm_prctile_hi(r);
            B_  = T.beta_median(r);
            Blo = T.beta_prctile_lo(r);
            Bhi = T.beta_prctile_hi(r);
            Sg  = T.sigma_median(r);

            if (Alo > 0) || (Ahi < 0), exc = 'YES'; else, exc = '-'; end

            fprintf('| %s -> %s | %.2f [%.2f, %.2f] | %.1f [%.1f, %.1f] | %.2f [%.2f, %.2f] | %.2f | %s |\n', ...
                labels{prev}, labels{curr}, A, Alo, Ahi, F, Flo, Fhi, B_, Blo, Bhi, Sg, exc);
        end
    end
end

%% Marginal summaries: pooled across the 9 cells per manipulation
fprintf('\n### Marginal summaries (median across the 9 cells per manipulation)\n\n');

for m = 1:2
    name = manip_labels{m};
    if m == 1, mtitle = 'Contrast'; else, mtitle = 'Precision'; end
    rows = T.cond_manipulation == name;
    fprintf('- **%s**: median A = %.2f deg (across-cell range %.2f-%.2f); median FWHM = %.1f deg (range %.1f-%.1f); median beta = %.2f deg.\n', ...
        mtitle, median(T.alpha_median(rows)), min(T.alpha_median(rows)), max(T.alpha_median(rows)), ...
        median(unbinnedWtoFwhm(T.w_median(rows))), ...
        min(unbinnedWtoFwhm(T.w_median(rows))), max(unbinnedWtoFwhm(T.w_median(rows))), ...
        median(T.beta_median(rows)));
end

fprintf('\n');
exit;
