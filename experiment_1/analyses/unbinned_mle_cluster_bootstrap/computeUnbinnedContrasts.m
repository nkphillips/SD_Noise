function ct = computeUnbinnedContrasts(results, contrast_specs)
% computeUnbinnedContrasts  Bootstrap CIs on linear-combination contrasts of the
% per-cell bootstrap parameters from runUnbinnedMLEClusterBootstrap.
%
% Inputs:
%   results        - results struct from runUnbinnedMLEClusterBootstrap. Must
%                    contain params_boot, jackknife.params, overlay.params_point,
%                    admitted, ci_prctile, B, n_subj.
%   contrast_specs - struct array with fields:
%                      .name       (char/string label)
%                      .weights    (1 x 18 double; linear combination over cells)
%                      .param_idx  (1=A, 2=w, 3=sigma, 4=beta)
%
% Output:
%   ct - table with one row per contrast: name, param_idx, estimate,
%        bca_lo, bca_hi, pc_lo, pc_hi, excludes_zero,
%        p_bca, p_bca_label, p_bonferroni, p_bonferroni_label,
%        p_holm, p_holm_label, p_fdr_bh, p_fdr_bh_label,
%        n_admit, z0, acc, ci_method.
%
% p_bca:           BCa-corrected two-sided bootstrap p-value (uncorrected for
%                  multiple comparisons).
% p_bonferroni:    p_bca * (number of contrasts in this table), capped at 1.
%                  Family = all rows in this contrast table. Conservative.
% p_holm:          Holm step-down adjustment, same family. Uniformly more
%                  powerful than Bonferroni and controls family-wise error.
% p_fdr_bh:        Benjamini-Hochberg FDR adjustment, same family. Controls
%                  expected false-discovery proportion rather than family-wise
%                  error; the modern default for exploratory work.
%
% Each *_label column is rounded to 1e-3 and decorated with stars:
%   ***  p < 0.001
%   **   p < 0.01
%   *    p < 0.05
%   (no star)  p >= 0.05
%
% IMPORTANT: the corrections use the FULL TABLE as the family. If you want to
% correct within a smaller family (e.g., only main effects of prev), subset the
% table first and recompute the corrections manually using the helpers below
% (or call computeUnbinnedContrasts on a subset of contrast_specs).
%
% A replicate contributes to a contrast only if all cells with non-zero weight
% in that contrast were admitted for that parameter.

ci_method = 'bca';
if isfield(results, 'ci_method') && ~isempty(results.ci_method)
    ci_method = lower(char(results.ci_method));
end
use_bca = strcmp(ci_method, 'bca');

if use_bca && (~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params'))
    error('computeUnbinnedContrasts:noJackknife', ...
        'BCa requires results.jackknife.params; rerun bootstrap with compute_jackknife = true.');
end

B = results.B;
n_subj = results.n_subj;
ci_pct = results.ci_prctile;
alpha = 1 - (ci_pct(2) - ci_pct(1)) / 100;
n_c = numel(contrast_specs);

name = strings(n_c, 1);
param_idx_col = zeros(n_c, 1);
estimate = nan(n_c, 1);
bca_lo = nan(n_c, 1); bca_hi = nan(n_c, 1);
pc_lo  = nan(n_c, 1); pc_hi  = nan(n_c, 1);
excludes_0 = false(n_c, 1);
p_bca_col = nan(n_c, 1);
p_bca_label_col = strings(n_c, 1);
n_admit_replicates = zeros(n_c, 1);
z0_col = nan(n_c, 1); acc_col = nan(n_c, 1);
ci_method_col = repmat(string(ci_method), n_c, 1);

for ic = 1:n_c
    spec = contrast_specs(ic);
    w_vec = reshape(spec.weights, 18, 1);    % 18 x 1 column
    k = spec.param_idx;
    invo = abs(w_vec) > sqrt(eps);           % logical 18 x 1

    name(ic) = string(spec.name);
    param_idx_col(ic) = k;

    % --- Bootstrap distribution: B x 1 contrast values ---
    boot_slice = reshape(results.params_boot(:, :, k), [B, 18]);  % B x 18
    % Replicate is good if every involved cell is admitted for this parameter.
    if isfield(results, 'param_admitted') && ~isempty(results.param_admitted)
        adm_slice = reshape(results.param_admitted(:, :, k), [B, 18]);
        adm_for_contrast = all(adm_slice(:, invo), 2);
    else
        adm_for_contrast = all(results.admitted(:, invo), 2);      % legacy saved result fallback
    end
    finite_for_contrast = all(isfinite(boot_slice(:, invo)), 2);   % B x 1 logical
    good = adm_for_contrast & finite_for_contrast;
    boot_v = nan(B, 1);
    boot_v(good) = boot_slice(good, :) * w_vec;                    % B x 18 * 18 x 1 = B x 1
    n_admit_replicates(ic) = sum(good);

    % --- Jackknife distribution: n_subj x 1 (BCa only) ---
    jk_v = nan(n_subj, 1);
    if use_bca
        jk_slice = reshape(results.jackknife.params(:, :, k), [n_subj, 18]);
        finite_jk = all(isfinite(jk_slice(:, invo)), 2);
        if isfield(results.jackknife, 'at_bound') && ~isempty(results.jackknife.at_bound) && ...
                isfield(results, 'discard_at_bound') && results.discard_at_bound
            jk_bound_slice = reshape(results.jackknife.at_bound(:, :, k), [n_subj, 18]);
            finite_jk = finite_jk & ~any(jk_bound_slice(:, invo), 2);
        end
        jk_v(finite_jk) = jk_slice(finite_jk, :) * w_vec;
    end

    % --- Point estimate from pooled fit ---
    pp_row = results.overlay.params_point(:, k);   % 18 x 1
    if all(isfinite(pp_row(invo)))
        estimate(ic) = pp_row' * w_vec;
    end

    bv = boot_v(isfinite(boot_v));
    if numel(bv) < 10
        continue
    end

    pc_lo(ic) = prctile(bv, ci_pct(1));
    pc_hi(ic) = prctile(bv, ci_pct(2));

    if use_bca
        [a1, a2, z0, a_acc] = bcaQuantilesLocal(boot_v, jk_v, estimate(ic), alpha);
        [bca_lo(ic), bca_hi(ic)] = orderedPrctileCI(bv, 100 * a1, 100 * a2);
        z0_col(ic) = z0;
        acc_col(ic) = a_acc;
        p_bca_col(ic) = bcaTwoSidedPValue(boot_v, z0, a_acc);
    else
        bca_lo(ic) = pc_lo(ic);
        bca_hi(ic) = pc_hi(ic);
        p_bca_col(ic) = percentileTwoSidedPValue(boot_v);
    end

    excludes_0(ic) = (bca_lo(ic) > 0) || (bca_hi(ic) < 0);

    p_bca_label_col(ic) = formatPLabel(p_bca_col(ic));
end

% --- Multiple-comparisons corrections, family = all rows in this table ---
p_bonf       = bonferroniAdjust(p_bca_col);
p_holm       = holmAdjust(p_bca_col);
p_bh         = benjaminiHochbergAdjust(p_bca_col);
p_bonf_label = arrayfun(@(x) formatPLabel(x), p_bonf);
p_holm_label = arrayfun(@(x) formatPLabel(x), p_holm);
p_bh_label   = arrayfun(@(x) formatPLabel(x), p_bh);

ct = table(name, param_idx_col, estimate, bca_lo, bca_hi, pc_lo, pc_hi, ...
    excludes_0, p_bca_col, p_bca_label_col, ...
    p_bonf, p_bonf_label, p_holm, p_holm_label, p_bh, p_bh_label, ...
    n_admit_replicates, z0_col, acc_col, ci_method_col, ...
    'VariableNames', {'name', 'param_idx', 'estimate', ...
    'bca_lo', 'bca_hi', 'pc_lo', 'pc_hi', ...
    'excludes_zero', 'p_bca', 'p_bca_label', ...
    'p_bonferroni', 'p_bonferroni_label', ...
    'p_holm', 'p_holm_label', ...
    'p_fdr_bh', 'p_fdr_bh_label', ...
    'n_admit', 'z0', 'acc', 'ci_method'});
end

function p_adj = bonferroniAdjust(p_raw)
% Bonferroni: p_adj = min(k * p_raw, 1).
k = numel(p_raw);
p_adj = min(k .* p_raw, 1);
p_adj(isnan(p_raw)) = NaN;
end

function p_adj = holmAdjust(p_raw)
% Holm step-down: sort ascending; multiplier (k - rank + 1) for rank-th smallest;
% enforce monotonicity (each adjusted p >= previous in sort order); cap at 1.
n = numel(p_raw);
p_adj = nan(n, 1);
valid = ~isnan(p_raw);
pv = p_raw(valid);
if isempty(pv)
    return
end
[pv_sorted, sort_idx] = sort(pv, 'ascend');
k = numel(pv_sorted);
multipliers = (k - (1:k) + 1)';   % column
p_step = pv_sorted .* multipliers;
% Enforce monotonicity along sort order
p_step = cummax(p_step);
p_step = min(p_step, 1);
p_resorted = nan(k, 1);
p_resorted(sort_idx) = p_step;
p_adj(valid) = p_resorted;
end

function p_adj = benjaminiHochbergAdjust(p_raw)
% Benjamini-Hochberg FDR: sort ascending; for rank i in 1..k,
%   q_(i) = min over j >= i of (p_(j) * k / j); cap at 1.
n = numel(p_raw);
p_adj = nan(n, 1);
valid = ~isnan(p_raw);
pv = p_raw(valid);
if isempty(pv)
    return
end
[pv_sorted, sort_idx] = sort(pv, 'ascend');
k = numel(pv_sorted);
ranks = (1:k)';
q_step = pv_sorted .* (k ./ ranks);
% Enforce monotonicity from largest rank toward smallest
q_step = flipud(cummin(flipud(q_step)));
q_step = min(q_step, 1);
q_resorted = nan(k, 1);
q_resorted(sort_idx) = q_step;
p_adj(valid) = q_resorted;
end

function p = bcaTwoSidedPValue(boot_v, z0, acc)
% Two-sided BCa p-value for H0: contrast == 0.
%
% F_BCa(v) = nominal alpha at which the BCa-adjusted alpha-quantile of the
%            bootstrap distribution equals v. Inverting:
%   q_v = empirical fraction of boot replicates < v
%   z_v = norminv(q_v)
%   s = z_v - z0
%   F_BCa(v) = normcdf(s/(1 + acc*s) - z0)
% Two-sided p = 2 * min(F_BCa(0), 1 - F_BCa(0)).
%
% Falls back to percentile-based two-sided p when BCa is numerically unstable.

bv = boot_v(isfinite(boot_v));
n = numel(bv);
if n < 10 || ~isfinite(z0) || ~isfinite(acc)
    p = NaN;
    return
end
q = mean(bv < 0);
% clamp away from {0, 1} to keep norminv finite
q = min(max(q, 1 / (n + 1)), 1 - 1 / (n + 1));
p_perc = 2 * min(q, 1 - q);

z_q = norminv(q);
s = z_q - z0;
denom = 1 + acc * s;
if denom <= 0 || ~isfinite(denom)
    p = p_perc;
    return
end
F_bca_0 = normcdf(s / denom - z0);
if ~isfinite(F_bca_0)
    p = p_perc;
    return
end
p = 2 * min(F_bca_0, 1 - F_bca_0);
p = min(max(p, 0), 1);
end

function p = percentileTwoSidedPValue(boot_v)
bv = boot_v(isfinite(boot_v));
n = numel(bv);
if n < 10
    p = NaN;
    return
end
q = mean(bv < 0);
q = min(max(q, 1 / (n + 1)), 1 - 1 / (n + 1));
p = 2 * min(q, 1 - q);
p = min(max(p, 0), 1);
end

function s = formatPLabel(p)
% Round p to 1e-3 and append significance stars.

if isnan(p)
    s = "NA";
    return
end
if p < 0.001
    stars = "***";
    body  = "<0.001";
elseif p < 0.01
    stars = "**";
    body  = string(sprintf('%.3f', round(p, 3)));
elseif p < 0.05
    stars = "*";
    body  = string(sprintf('%.3f', round(p, 3)));
else
    stars = "";
    body  = string(sprintf('%.3f', round(p, 3)));
end
s = body + stars;
end

function [lo, hi] = orderedPrctileCI(values, prc_lo, prc_hi)
q1 = prctile(values, prc_lo);
q2 = prctile(values, prc_hi);
lo = min(q1, q2);
hi = max(q1, q2);
end

function [a1, a2, z0, a_acc] = bcaQuantilesLocal(boot_vals, jk_vals, point_est, alpha)
boot_vals = boot_vals(isfinite(boot_vals));
jk_vals = jk_vals(isfinite(jk_vals));
a1 = alpha / 2; a2 = 1 - alpha / 2; z0 = NaN; a_acc = NaN;
if numel(boot_vals) < 2 || numel(jk_vals) < 2 || ~isfinite(point_est)
    return
end
p_hat = mean(boot_vals < point_est);
p_hat = min(max(p_hat, 1 / (numel(boot_vals) + 1)), 1 - 1 / (numel(boot_vals) + 1));
z0 = norminv(p_hat);
jk_mean = mean(jk_vals);
num_a = sum((jk_mean - jk_vals).^3);
den_a = 6 * (sum((jk_mean - jk_vals).^2)).^(3/2);
if den_a == 0
    a_acc = 0;
else
    a_acc = num_a / den_a;
end
z_lo = norminv(alpha / 2);
z_hi = norminv(1 - alpha / 2);
a1 = normcdf(z0 + (z0 + z_lo) / (1 - a_acc * (z0 + z_lo)));
a2 = normcdf(z0 + (z0 + z_hi) / (1 - a_acc * (z0 + z_hi)));
if ~isfinite(a1) || a1 <= 0, a1 = alpha / 2; end
if ~isfinite(a2) || a2 >= 1, a2 = 1 - alpha / 2; end
end
