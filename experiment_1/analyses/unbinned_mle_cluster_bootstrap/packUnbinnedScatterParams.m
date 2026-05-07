function pack = packUnbinnedScatterParams(summary_table)
% packUnbinnedScatterParams  Build 3×3×2 grids from results.summary_table for scatter plots.
%
% pack.amp, pack.fwhm — point estimates (median bootstrap α, FWHM from w_median)
% pack.amp_lo/hi, pack.fwhm_lo/hi — bootstrap percentile bounds on α and on FWHM(w).

    pack.amp = nan(3, 3, 2);
    pack.fwhm = nan(3, 3, 2);
    pack.amp_lo = nan(3, 3, 2);
    pack.amp_hi = nan(3, 3, 2);
    pack.fwhm_lo = nan(3, 3, 2);
    pack.fwhm_hi = nan(3, 3, 2);

    for r = 1:height(summary_table)
        prev = summary_table.cond_prev(r);
        curr = summary_table.cond_curr(r);
        if summary_table.cond_manipulation(r) == 'contrast'
            m = 1;
        else
            m = 2;
        end

        pack.amp(prev, curr, m) = summary_table.alpha_median(r);
        al = summary_table.alpha_prctile_lo(r);
        ah = summary_table.alpha_prctile_hi(r);
        pack.amp_lo(prev, curr, m) = min(al, ah);
        pack.amp_hi(prev, curr, m) = max(al, ah);

        w0 = summary_table.w_median(r);
        pack.fwhm(prev, curr, m) = unbinnedWtoFwhm(w0);

        wl = summary_table.w_prctile_lo(r);
        wh = summary_table.w_prctile_hi(r);
        % FWHM = 2*sqrt(log(2))/w (decreasing in w under S&S param); min/max ensures lo <= hi
        fw1 = unbinnedWtoFwhm(wl);
        fw2 = unbinnedWtoFwhm(wh);
        pack.fwhm_lo(prev, curr, m) = min(fw1, fw2);
        pack.fwhm_hi(prev, curr, m) = max(fw1, fw2);
    end

end
