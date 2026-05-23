function pack = packSerialDependenceScatterParams(summary_table)
% packSerialDependenceScatterParams  Build 3×3×2 grids from results.summary_table for scatter plots.
%
% pack.amp, pack.fwhm — point estimates from the pooled full-data fit.
% pack.amp_lo/hi, pack.fwhm_lo/hi — active bootstrap CI bounds from summary_table.

    pack.amp = nan(3, 3, 2);
    pack.fwhm = nan(3, 3, 2);
    pack.amp_lo = nan(3, 3, 2);
    pack.amp_hi = nan(3, 3, 2);
    pack.fwhm_lo = nan(3, 3, 2);
    pack.fwhm_hi = nan(3, 3, 2);

    for r = 1:height(summary_table)
        prev = summary_table.cond_prev(r);
        curr = summary_table.cond_curr(r);
        if strcmpi(char(summary_table.cond_manipulation(r)), 'contrast')
            m = 1;
        else
            m = 2;
        end

        pack.amp(prev, curr, m) = summary_table.A_point(r);
        al = summary_table.A_ci_lo(r);
        ah = summary_table.A_ci_hi(r);
        pack.amp_lo(prev, curr, m) = min(al, ah);
        pack.amp_hi(prev, curr, m) = max(al, ah);

        pack.fwhm(prev, curr, m) = summary_table.fwhm_point(r);

        fl = summary_table.fwhm_ci_lo(r);
        fh = summary_table.fwhm_ci_hi(r);
        pack.fwhm_lo(prev, curr, m) = min(fl, fh);
        pack.fwhm_hi(prev, curr, m) = max(fl, fh);
    end

end
