function overlay = computePooledFitsAndOverlay(tbl_trials, delta_centers, window_width, guess_rate)
% Per-condition pooled MLE on full tbl_trials (S&S parameterization, 25% guess rate),
% fit metrics vs predicted P(CW), binned μ MLE for figures.

    if nargin < 4 || isempty(guess_rate)
        guess_rate = 0.25;
    end

    num_conds = 18;
    overlay.params_point = nan(num_conds, 4);
    overlay.r2 = nan(num_conds, 1);           % legacy: same as r2_efron (Efron R² on binary trials)
    overlay.r2_efron = nan(num_conds, 1);
    overlay.r2_tjur = nan(num_conds, 1);
    overlay.r2_mcfadden = nan(num_conds, 1);
    overlay.r2_delta_bins = nan(num_conds, 1);   % Efron on Δθ-bin means (weighted)—aligned with figure
    overlay.mu_bin_mle = cell(num_conds, 1);
    overlay.guess_rate = guess_rate;

    for c = 1:num_conds
        [m, prev, curr] = conditionSubscriptsFromIndex(c);
        cm = tbl_trials.cond_manipulation;
        man = ones(height(tbl_trials), 1);
        man(cm == 'precision') = 2;
        mask = man == m & tbl_trials.cond_prev == prev & tbl_trials.cond_curr == curr;

        dt = tbl_trials.delta_theta(mask);
        xp = tbl_trials.x_probe(mask);
        rv = tbl_trials.response(mask);

        pfit = fitConditionMLE(dt, xp, rv, struct('guess_rate', guess_rate));
        overlay.params_point(c, :) = pfit;

        if all(isnan(pfit))
            continue
        end

        mu_t = dogIsolated(dt, pfit(1), pfit(2)) + pfit(4);
        p_psy = normcdf(xp, mu_t, pfit(3));
        p_hat = (1 - guess_rate) .* p_psy + 0.5 .* guess_rate;
        met = binaryFitR2Metrics(rv(:), p_hat(:));
        overlay.r2_efron(c) = met.efron;
        overlay.r2_tjur(c) = met.tjur;
        overlay.r2_mcfadden(c) = met.mcfadden;
        overlay.r2(c) = met.efron;

        overlay.r2_delta_bins(c) = binnedEfronAlongDeltaTheta(dt, rv(:), p_hat(:), delta_centers, window_width);

        [dc, mu_b] = mleMuPerDeltaBin(dt, xp, rv, pfit(3), delta_centers, window_width, guess_rate);
        overlay.mu_bin_mle{c} = struct('delta_centers', dc(:), 'mu_deg', mu_b(:));
    end

end
