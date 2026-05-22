function w = unbinnedFwhmToW(fwhm)
% unbinnedFwhmToW  Convert S&S DoG-lobe FWHM in deg to width w (1/deg).

    fwhm = fwhm(:);
    w = nan(size(fwhm));
    ok = isfinite(fwhm) & fwhm > 0;
    if any(ok)
        coeff = unbinnedWtoFwhm(1);
        w(ok) = coeff ./ fwhm(ok);
    end
end
