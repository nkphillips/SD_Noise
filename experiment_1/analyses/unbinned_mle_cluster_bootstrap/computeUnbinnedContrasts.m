function ct = computeUnbinnedContrasts(results, contrast_specs)
% computeUnbinnedContrasts  BCa CIs on linear-combination contrasts of the
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
%        bca_lo, bca_hi, pc_lo, pc_hi, excludes_zero, n_admit, z0, acc.
%
% A replicate contributes to a contrast only if all cells with non-zero weight
% in that contrast were admitted in that replicate (admission filter inherited
% from the BCa pipeline -- discard_at_bound and discard_failed_fits).

    if ~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params')
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
    n_admit_replicates = zeros(n_c, 1);
    z0_col = nan(n_c, 1); acc_col = nan(n_c, 1);

    for ic = 1:n_c
        spec = contrast_specs(ic);
        w_vec = reshape(spec.weights, 18, 1);    % 18 x 1 column
        k = spec.param_idx;
        invo = abs(w_vec) > sqrt(eps);           % logical 18 x 1

        name(ic) = string(spec.name);
        param_idx_col(ic) = k;

        % --- Bootstrap distribution: B x 1 contrast values ---
        boot_slice = reshape(results.params_boot(:, :, k), [B, 18]);  % B x 18
        % admitted cells in each replicate involved in the contrast
        % (replicate is "good" for this contrast iff every involved cell is admitted)
        adm_for_contrast = all(results.admitted(:, invo), 2);          % B x 1 logical
        finite_for_contrast = all(isfinite(boot_slice(:, invo)), 2);   % B x 1 logical
        good = adm_for_contrast & finite_for_contrast;
        boot_v = nan(B, 1);
        boot_v(good) = boot_slice(good, :) * w_vec;                    % B x 18 * 18 x 1 = B x 1
        n_admit_replicates(ic) = sum(good);

        % --- Jackknife distribution: n_subj x 1 ---
        jk_slice = reshape(results.jackknife.params(:, :, k), [n_subj, 18]);
        finite_jk = all(isfinite(jk_slice(:, invo)), 2);
        jk_v = nan(n_subj, 1);
        jk_v(finite_jk) = jk_slice(finite_jk, :) * w_vec;

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

        [a1, a2, z0, a_acc] = bcaQuantilesLocal(boot_v, jk_v, estimate(ic), alpha);
        bca_lo(ic) = prctile(bv, 100 * a1);
        bca_hi(ic) = prctile(bv, 100 * a2);
        z0_col(ic) = z0;
        acc_col(ic) = a_acc;

        excludes_0(ic) = (bca_lo(ic) > 0) || (bca_hi(ic) < 0);
    end

    ct = table(name, param_idx_col, estimate, bca_lo, bca_hi, pc_lo, pc_hi, ...
               excludes_0, n_admit_replicates, z0_col, acc_col, ...
               'VariableNames', {'name', 'param_idx', 'estimate', ...
                                 'bca_lo', 'bca_hi', 'pc_lo', 'pc_hi', ...
                                 'excludes_zero', 'n_admit', 'z0', 'acc'});
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
